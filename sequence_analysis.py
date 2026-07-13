"""
Dengue Genome Aggregator - Single-Sequence Analysis Pipeline

The read-only counterpart to the batch ingestion pipeline (main.py /
annotate_variants.py): given one pasted/uploaded DENV sequence (complete
genome or partial fragment), determine which serotype it is, call variants
against the matching RefSeq reference, cross-reference those variants
against everything already known in dengue_genomes.db (population
frequency, geography, epitope overlap), and place it against its nearest
neighbors in the existing 12k-genome collection.

Nothing here writes to dengue_genomes.db -- every connection is opened with
PRAGMA query_only, so any accidental write raises instead of silently
succeeding. See fingerprints.py for the nearest-neighbor precomputation this
depends on (run build_fingerprints.py first).
"""

import logging
import re
import sqlite3
from pathlib import Path
from typing import Optional

import numpy as np
from rapidfuzz.distance import Levenshtein

import fingerprints
import references
from variant_annotator import VariantAnnotator

logger = logging.getLogger(__name__)

DB_PATH = Path(__file__).parent / "dengue_genomes.db"

IUPAC_NUCLEOTIDES = set("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

KMER_SIZE = 21
# Even a correct serotype call rarely gets anywhere near 100% k-mer
# containment: real genomes diverge a few percent from the RefSeq reference
# *within* their own serotype, and at k=21 a single mismatch kills 21
# k-mers, so e.g. 7% intra-serotype divergence knocks containment down to
# roughly (1-0.07)**21 =~ 21% even for the right call. Serotypes diverge
# ~30-35% from each other, so the signal that actually separates a correct
# call isn't the absolute score -- it's how far ahead of the runner-up
# serotype the winner is.
MIN_SIGNAL_PCT = 2.0      # winner below this: essentially no k-mer signal at all
MIN_MARGIN_RATIO = 3.0    # winner must beat the runner-up by at least this factor
BASE_CODE = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

_catalog_cache: dict = {}
_epitope_cache: Optional[list] = None
_reference_kmer_cache: dict = {}


class SequenceInputError(ValueError):
    """Raised for input that can't be analyzed at all (empty, not nucleotide data, etc.)."""


def _get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA query_only = TRUE")
    conn.row_factory = sqlite3.Row
    return conn


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _clean_input(raw_text: str) -> tuple[str, list[str]]:
    """
    Strip FASTA header(s) and whitespace from pasted or uploaded text.
    Only the first FASTA record (or the whole input, if headerless) is
    analyzed; additional records produce a warning rather than being silently
    dropped.

    Deliberately does *not* parse every record just to discard all but the
    first: a user pointing this at a real multi-genome batch export (e.g.
    one of the Nextstrain per-serotype FASTAs in dengue_data/, which run to
    100MB+ for thousands of genomes) would otherwise pay to fully tokenize
    the whole file for nothing. `str.find`/`str.count` are C-level scans, so
    slicing out just the first record and counting the rest stays fast
    regardless of total input size.
    """
    warnings = []
    text = raw_text.strip()
    if not text:
        raise SequenceInputError("No sequence provided.")

    if text.startswith('>'):
        header_end = text.find('\n')
        header = text[1:header_end].strip() if header_end != -1 else text[1:].strip()

        if header_end == -1:
            body, rest = '', ''
        else:
            next_header = text.find('\n>', header_end)
            body = text[header_end + 1:] if next_header == -1 else text[header_end + 1:next_header]
            rest = '' if next_header == -1 else text[next_header + 1:]

        record_count = 1 + rest.count('\n>') + (1 if rest.startswith('>') else 0)
        seq = body
    else:
        header = None
        record_count = 1
        seq = text

    seq = re.sub(r'[^A-Za-z]', '', seq).upper()
    if not seq:
        raise SequenceInputError("Could not find any sequence data in the input.")

    if record_count > 1:
        warnings.append(
            f"Input contained {record_count} FASTA records; only the first "
            f"(\"{(header or '')[:60]}\") was analyzed."
        )

    return seq, warnings


def _validate_and_clean(seq: str) -> tuple[str, list[str]]:
    """Check alphabet/length; replace stray non-IUPAC characters with N."""
    warnings = []
    if len(seq) < 30:
        raise SequenceInputError(
            f"Sequence is too short to analyze ({len(seq)} bp; need at least 30 bp)."
        )

    invalid_chars = set(seq) - IUPAC_NUCLEOTIDES
    if invalid_chars:
        invalid_fraction = sum(seq.count(c) for c in invalid_chars) / len(seq)
        if invalid_fraction > 0.05:
            raise SequenceInputError(
                f"{invalid_fraction*100:.1f}% of characters aren't valid nucleotide codes "
                f"(found: {''.join(sorted(invalid_chars))[:20]}). Is this a protein sequence, "
                f"or does the pasted text contain non-sequence content?"
            )
        warnings.append(
            f"{invalid_fraction*100:.1f}% of characters were non-IUPAC and were treated as N."
        )
        seq = ''.join(c if c in IUPAC_NUCLEOTIDES else 'N' for c in seq)

    n_count = seq.count('N')
    if len(seq) and n_count / len(seq) > 0.10:
        warnings.append(f"High ambiguous-base content: {n_count/len(seq)*100:.1f}% N.")

    return seq, warnings


def _kmer_set(seq: str, k: int) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1) if set(seq[i:i + k]) <= set('ACGT')}


def _reference_kmers(serotype: int, k: int) -> set:
    cache_key = (serotype, k)
    if cache_key not in _reference_kmer_cache:
        ref_seq = references.get_reference(serotype)['sequence']
        _reference_kmer_cache[cache_key] = _kmer_set(ref_seq, k)
    return _reference_kmer_cache[cache_key]


_reference_kmer_position_cache: dict = {}


def _reference_kmer_positions(serotype: int, k: int) -> dict:
    """kmer -> list of 0-indexed positions in the serotype reference (used for seeding)."""
    cache_key = (serotype, k)
    if cache_key not in _reference_kmer_position_cache:
        ref_seq = references.get_reference(serotype)['sequence']
        positions: dict = {}
        for i in range(len(ref_seq) - k + 1):
            kmer = ref_seq[i:i + k]
            if set(kmer) <= set('ACGT'):
                positions.setdefault(kmer, []).append(i)
        _reference_kmer_position_cache[cache_key] = positions
    return _reference_kmer_position_cache[cache_key]


def _estimate_offset(serotype: int, query_seq: str, k: int = KMER_SIZE) -> Optional[int]:
    """
    Estimate the 0-indexed reference position the query starts at, via exact
    k-mer seed matching (mode of per-seed offset votes). Each surviving exact
    match gives a deterministic, driftless offset when there's no indel
    between reference and query up to that point, so the mode is normally
    exactly right, not just approximately. Returns None if too few seeds
    survive to localize (very short or very divergent query).
    """
    if len(query_seq) < k:
        return None

    positions_index = _reference_kmer_positions(serotype, k)
    step = max(1, k // 3)
    votes = []
    for i in range(0, len(query_seq) - k + 1, step):
        kmer = query_seq[i:i + k]
        for pos in positions_index.get(kmer, ()):
            votes.append(pos - i)

    if not votes:
        return None

    from collections import Counter
    estimated_offset, _ = Counter(votes).most_common(1)[0]
    return estimated_offset


def align_query(serotype: int, query_seq: str):
    """
    Localize the query against the serotype reference and align it, choosing
    a reference window deliberately sized to avoid a specific failure mode of
    plain global edit distance: aligning a short fragment against the *full*
    (much longer) reference degenerates into scattering the query across
    non-contiguous matching bases anywhere in the genome, because with only 4
    letters in the alphabet, exploiting a coincidental match is essentially
    always "free" whenever there's spare deletion budget to shuffle around
    (confirmed empirically -- a real 1.5kb E-gene fragment aligned against
    the full ~10.7kb reference came back with 0 SNPs and 1074 scattered
    single-base indels covering 55% of the genome, instead of the ~15 real
    SNPs confined to the true E-gene window; even a generously-padded 300bp
    window still leaked ~150bp of this artifact into the reported coverage
    span). The fix: window the reference as tightly as possible.

    Strategy: seed-localize the query's start offset, then try a zero-slack
    window (exactly query-length, eliminating the free-deletion budget the
    scattering exploit needs) first. That's exact whenever there's no net
    indel between reference and query before this point. If the resulting
    edit rate is implausibly high (>25%, suggesting a real indel shifted the
    reading frame), retry with a small amount of padding to give the aligner
    room to place that indel, and keep whichever attempt has the lower edit
    rate.

    Returns (window, ops, offset_found) where window is a (start, end)
    0-indexed reference slice and ops is the rapidfuzz editops result for
    that window vs. query_seq. offset_found is False if seeding failed
    entirely, in which case window falls back to the full reference.
    """
    ref = references.get_reference(serotype)
    ref_len = ref['length']
    ref_seq = ref['sequence']

    offset = _estimate_offset(serotype, query_seq)
    if offset is None:
        window = (0, ref_len)
        ops = Levenshtein.editops(ref_seq, query_seq)
        return window, ops, False

    def _try(pad: int):
        start = max(0, offset - pad)
        end = min(ref_len, offset + len(query_seq) + pad)
        ops = Levenshtein.editops(ref_seq[start:end], query_seq)
        edit_rate = len(ops) / len(query_seq) if query_seq else 1.0
        return (start, end), ops, edit_rate

    window, ops, edit_rate = _try(pad=0)
    if edit_rate > 0.25:
        window2, ops2, edit_rate2 = _try(pad=30)
        if edit_rate2 < edit_rate:
            window, ops = window2, ops2

    return window, ops, True


def _kmer_containment(query_kmers: set, ref_kmers: set) -> float:
    if not query_kmers:
        return 0.0
    return 100.0 * len(query_kmers & ref_kmers) / len(query_kmers)


def call_serotype_and_orientation(seq: str) -> dict:
    """
    Determine serotype + strand orientation via exact k-mer containment
    against the 4 RefSeq references, tried in both the given orientation and
    its reverse complement (a pasted/uploaded sequence's strand isn't known
    up front). DENV serotypes diverge ~30-35% at the nucleotide level, so
    containment scores separate cleanly for real data.

    Falls back to a shorter k for fragments smaller than one k-mer.
    """
    k = KMER_SIZE if len(seq) >= KMER_SIZE else max(8, len(seq) // 2 or 1)
    seq_rc = reverse_complement(seq)

    fwd_kmers = _kmer_set(seq, k)
    rc_kmers = _kmer_set(seq_rc, k)

    scores = {}
    for serotype in (1, 2, 3, 4):
        ref_kmers = _reference_kmers(serotype, k)
        scores[(serotype, 'forward')] = _kmer_containment(fwd_kmers, ref_kmers)
        scores[(serotype, 'reverse_complement')] = _kmer_containment(rc_kmers, ref_kmers)

    best_key = max(scores, key=scores.get)
    best_serotype, best_orientation = best_key
    best_score = scores[best_key]

    # Report all 4 serotypes' scores in whichever orientation won, so the UI
    # can show the evidence behind the call, not just the winner.
    scores_by_serotype = {
        st: round(scores[(st, best_orientation)], 1) for st in (1, 2, 3, 4)
    }
    ranked = sorted(scores_by_serotype.values(), reverse=True)
    runner_up = ranked[1] if len(ranked) > 1 else 0.0
    margin_ratio = (ranked[0] / runner_up) if runner_up > 0 else float('inf')
    no_signal = ranked[0] < MIN_SIGNAL_PCT
    ambiguous = no_signal or margin_ratio < MIN_MARGIN_RATIO

    return {
        'serotype': best_serotype,
        'orientation': best_orientation,
        'confidence_pct': round(best_score, 1),
        'scores_by_serotype': scores_by_serotype,
        'margin_ratio': round(margin_ratio, 1) if margin_ratio != float('inf') else None,
        'no_signal': no_signal,
        'ambiguous': ambiguous,
    }


def _covered_span(ref_len: int, ops) -> tuple[int, int]:
    """
    Derive the reference region actually spanned by the query from the
    editops list: a run of 'delete' ops anchored at reference position 0
    means the query doesn't start until later (a fragment, not a genome with
    a real 5' deletion -- see sequence_analysis's module docstring / the
    design plan for why that ambiguity is accepted). Same logic in reverse
    for the 3' end. Returns 0-indexed [start, end).
    """
    lead = 0
    while lead < len(ops) and ops[lead].tag == 'delete' and ops[lead].src_pos == lead:
        lead += 1

    trail = 0
    j = len(ops) - 1
    while j >= 0 and ops[j].tag == 'delete' and ops[j].src_pos == ref_len - 1 - trail:
        trail += 1
        j -= 1

    return lead, ref_len - trail


def _load_catalog(serotype: int) -> dict:
    if serotype not in _catalog_cache:
        conn = _get_connection()
        try:
            rows = conn.execute(
                "SELECT * FROM variant_catalog WHERE serotype = ?", (serotype,)
            ).fetchall()
        finally:
            conn.close()
        _catalog_cache[serotype] = {
            (r['position'], r['ref_allele'], r['alt_allele']): dict(r) for r in rows
        }
        logger.info(f"Cached variant_catalog for serotype {serotype}: {len(rows)} rows")
    return _catalog_cache[serotype]


def _load_epitopes() -> list:
    global _epitope_cache
    if _epitope_cache is None:
        conn = _get_connection()
        try:
            rows = conn.execute("SELECT * FROM epitopes").fetchall()
        finally:
            conn.close()
        _epitope_cache = [dict(r) for r in rows]
        logger.info(f"Cached epitopes: {len(_epitope_cache)} rows")
    return _epitope_cache


def analyze_sequence(raw_text: str) -> dict:
    """Run the full pipeline on pasted/uploaded text. Raises SequenceInputError
    for unusable input; any other exception indicates a real bug."""
    warnings: list[str] = []

    seq, parse_warnings = _clean_input(raw_text)
    warnings += parse_warnings

    seq, validate_warnings = _validate_and_clean(seq)
    warnings += validate_warnings

    length = len(seq)
    gc_count = seq.count('G') + seq.count('C')
    n_count = seq.count('N')

    serotype_call = call_serotype_and_orientation(seq)

    if serotype_call['no_signal']:
        warnings.append(
            f"Almost no k-mer matches to any DENV reference (best: DENV-{serotype_call['serotype']} at "
            f"{serotype_call['confidence_pct']}%) -- this likely isn't a dengue sequence, or is "
            "extremely degraded/short."
        )
    elif serotype_call['ambiguous']:
        warnings.append(
            f"Serotype call isn't decisively separated from the runner-up (DENV-{serotype_call['serotype']} "
            f"leads by only {serotype_call['margin_ratio']}x) -- consider a mixed infection, contamination, "
            "or a recombinant/highly divergent sequence."
        )

    serotype = serotype_call['serotype']
    query_seq = seq if serotype_call['orientation'] == 'forward' else reverse_complement(seq)

    ref = references.get_reference(serotype)
    ref_len = ref['length']

    window, ops, offset_found = align_query(serotype, query_seq)
    if not offset_found:
        warnings.append(
            "Could not localize this sequence against the reference via k-mer seeding "
            "(too short or too divergent) -- falling back to whole-genome alignment, which "
            "may misplace variants for a partial sequence."
        )

    annotator = VariantAnnotator(serotype)
    variants = annotator.call_variants(query_seq, ref_window=window)

    win_start, win_end = window
    local_start, local_end = _covered_span(win_end - win_start, ops)
    cov_start, cov_end = win_start + local_start, win_start + local_end
    coverage_pct = round(100 * (cov_end - cov_start) / ref_len, 1) if ref_len else 0.0

    if coverage_pct < 95:
        warnings.append(
            f"Query aligns to {coverage_pct}% of the DENV-{serotype} reference genome "
            f"(positions {cov_start + 1}-{cov_end} of {ref_len}) -- this is a partial sequence; "
            "variant/lineage calls below only reflect that region."
        )

    catalog = _load_catalog(serotype)
    epitope_rows = _load_epitopes()

    variant_dicts = []
    novel_count = 0
    epitope_overlap_count = 0
    impact_counts = {'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'MODIFIER': 0}
    type_counts = {'SNP': 0, 'INDEL': 0}

    query_vec = np.zeros(ref_len, dtype=np.uint8)

    for v in variants:
        vd = annotator.annotate_with_epitopes(v, epitope_rows)

        key = (v.position, v.ref_allele, v.alt_allele)
        known = catalog.get(key)
        if known:
            vd['known'] = True
            vd['occurrence_count'] = known['occurrence_count']
            vd['genome_count_serotype'] = known['genome_count_serotype']
            vd['frequency_pct'] = round((known['frequency'] or 0) * 100, 3)
            vd['first_observed_date'] = known['first_observed_date']
            vd['last_observed_date'] = known['last_observed_date']
            vd['countries'] = known['countries']
            vd['country_count'] = known['country_count']
            vd['lineages'] = known['lineages']
        else:
            vd['known'] = False
            novel_count += 1

        vd['low_confidence'] = v.variant_type == 'SNP' and v.alt_allele not in 'ACGT'

        if vd['epitope_hits']:
            epitope_overlap_count += 1

        impact_counts[v.impact] = impact_counts.get(v.impact, 0) + 1
        type_counts[v.variant_type] = type_counts.get(v.variant_type, 0) + 1

        if v.variant_type == 'SNP' and v.alt_allele in BASE_CODE:
            query_vec[v.position - 1] = BASE_CODE[v.alt_allele]

        variant_dicts.append(vd)

    variant_dicts.sort(key=lambda d: d['position'])

    covered_mask = np.zeros(ref_len, dtype=bool)
    covered_mask[cov_start:cov_end] = True

    try:
        neighbors = fingerprints.nearest_neighbors(serotype, query_vec, covered_mask, top_n=10)
    except FileNotFoundError as e:
        neighbors = []
        warnings.append(str(e))

    return {
        'input_summary': {
            'length': length,
            'gc_pct': round(100 * gc_count / length, 1) if length else 0.0,
            'n_count': n_count,
            'n_pct': round(100 * n_count / length, 1) if length else 0.0,
        },
        'serotype_call': serotype_call,
        'coverage': {
            'start': cov_start + 1,
            'end': cov_end,
            'pct': coverage_pct,
            'reference_length': ref_len,
        },
        'reference_info': {
            'serotype': serotype,
            'accession': ref['accession'],
            'length': ref_len,
            'proteins': ref['proteins'],
        },
        'variants': variant_dicts,
        'summary_counts': {
            'total': len(variant_dicts),
            'by_impact': impact_counts,
            'by_type': type_counts,
            'novel': novel_count,
            'epitope_overlap': epitope_overlap_count,
        },
        'nearest_neighbors': neighbors,
        'warnings': warnings,
    }
