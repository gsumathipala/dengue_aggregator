"""
Dengue Genome Aggregator - Epitope Loader

Loads the IEDB linear B-/T-cell epitope records fetched by fetch_iedb.py
(references/iedb_dengue_linear_epitopes.json) into the `epitopes` table, and
back-fills the epitope-overlap columns on `variant_catalog` (epitope_overlap,
epitope_hits) so the immune-escape-risk queries in the README actually
return rows. Run after fetch_iedb.py and (for the variant_catalog backfill
step) build_variant_catalog.py.

IEDB doesn't report positions in our reference genomes' numbering -- its
`curated_source_antigens` positions are amino-acid offsets into whichever
isolate's polyprotein/protein record was curated, which don't line up with
our NC_00147x coordinates. Instead, each epitope's peptide sequence is
located by exact substring search against the mature-protein amino acid
sequences translated from our own serotype reference genomes (references.py),
which is exactly the coordinate system variant_annotator.py's aa_change /
protein_position values already use -- so a match here is directly joinable
against `variants.protein` + `variants.protein_position`.

Usage:
    python load_epitopes.py --db dengue_genomes.db
"""

import argparse
import json
import logging
import re
import uuid
from collections import defaultdict
from pathlib import Path

from db_schema import DengueDB
from variant_annotator import GENETIC_CODE
import references

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

IEDB_JSON_PATH = Path(__file__).parent / "references" / "iedb_dengue_linear_epitopes.json"

# NCBI Taxonomy IDs for the four DENV serotypes, cross-checked in references.py.
TAXON_TO_SEROTYPE = {
    "NCBITaxon:11053": 1,
    "NCBITaxon:11060": 2,
    "NCBITaxon:11069": 3,
    "NCBITaxon:11070": 4,
}

MIN_PEPTIDE_LENGTH = 5  # shorter peptides match too many spurious positions to be useful


def translate(nt_seq: str) -> str:
    """Translate a nucleotide sequence to amino acids, one triplet at a time."""
    aa = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i:i + 3]
        aa.append(GENETIC_CODE.get(codon, 'X'))
    return ''.join(aa)


def build_protein_aa_sequences() -> dict:
    """{serotype: {protein_name: aa_sequence}} translated from the reference genomes."""
    result = {}
    for serotype in (1, 2, 3, 4):
        ref = references.get_reference(serotype)
        result[serotype] = {
            name: translate(ref['sequence'][coords['start']:coords['end']])
            for name, coords in ref['proteins'].items()
        }
    return result


def candidate_serotypes(antigens: list) -> list:
    """Serotypes suggested by an epitope's curated_source_antigens; empty -> search all."""
    serotypes = set()
    for antigen in antigens or []:
        iri = antigen.get('source_organism_iri')
        if iri in TAXON_TO_SEROTYPE:
            serotypes.add(TAXON_TO_SEROTYPE[iri])
    return sorted(serotypes)


def find_positions(peptide: str, protein_seqs: dict, serotypes: list) -> set:
    """
    Search `peptide` in the given serotypes' translated proteins.

    Returns a set of (protein, start_pos, end_pos), 1-indexed amino-acid
    coordinates within that protein -- matching the units of
    variants.protein_position (see variant_annotator.VariantAnnotator).
    """
    hits = set()
    for serotype in serotypes:
        for protein_name, aa_seq in protein_seqs[serotype].items():
            start = 0
            while True:
                idx = aa_seq.find(peptide, start)
                if idx == -1:
                    break
                hits.add((protein_name, idx + 1, idx + len(peptide)))
                start = idx + 1
    return hits


def classify_epitope_type(record: dict) -> str:
    has_antibody = bool(record.get('antibody_isotypes'))
    has_mhc = bool(record.get('mhc_allele_names'))
    if has_antibody and has_mhc:
        return 'B-cell/T-cell'
    if has_antibody:
        return 'B-cell'
    if has_mhc:
        return 'T-cell'
    return 'Unknown'


def summarize_immune_response(record: dict) -> str:
    quals = record.get('qualitative_measures') or []
    isotypes = record.get('antibody_isotypes') or []
    parts = sorted(set(quals)) + sorted(set(isotypes))
    return '; '.join(parts) if parts else None


def load_epitopes(db_path: str, json_path: Path = IEDB_JSON_PATH) -> dict:
    if not json_path.exists():
        raise FileNotFoundError(f"IEDB epitope file not found: {json_path}. Run fetch_iedb.py first.")

    with open(json_path, encoding='utf-8') as f:
        records = json.load(f)
    logger.info(f"Loaded {len(records)} IEDB epitope records from {json_path.name}")

    protein_seqs = build_protein_aa_sequences()
    all_serotypes = [1, 2, 3, 4]

    stats = {'total': len(records), 'too_short': 0, 'unmapped': 0, 'mapped_records': 0, 'rows_inserted': 0}
    insert_rows = []

    for record in records:
        peptide = (record.get('linear_sequence') or '').upper()
        if len(peptide) < MIN_PEPTIDE_LENGTH or not re.fullmatch(r'[A-Z]+', peptide):
            stats['too_short'] += 1
            continue

        antigens = record.get('curated_source_antigens') or []
        serotypes = candidate_serotypes(antigens) or all_serotypes

        hits = find_positions(peptide, protein_seqs, serotypes)
        if not hits and serotypes != all_serotypes:
            # Curated antigen's serotype didn't line up with our reference (isolate
            # divergence) -- fall back to searching every serotype before giving up.
            hits = find_positions(peptide, protein_seqs, all_serotypes)

        if not hits:
            stats['unmapped'] += 1
            continue

        stats['mapped_records'] += 1
        epitope_type = classify_epitope_type(record)
        immune_response = summarize_immune_response(record)

        for protein_name, start_pos, end_pos in hits:
            insert_rows.append((
                f"IEDB_{record['structure_id']}_{uuid.uuid4().hex[:8]}",
                protein_name, start_pos, end_pos, peptide,
                epitope_type, immune_response, 'IEDB',
            ))

    stats['rows_inserted'] = len(insert_rows)
    logger.info(f"Mapped {stats['mapped_records']}/{stats['total']} epitopes to a reference protein "
                f"position ({stats['unmapped']} unmapped, {stats['too_short']} too short/non-peptide) "
                f"-> {stats['rows_inserted']} epitope rows (an epitope can map to more than one position)")

    db = DengueDB(db_path)
    cur = db.conn.cursor()
    cur.execute("DELETE FROM epitopes WHERE source_database = 'IEDB'")
    cur.executemany("""
        INSERT INTO epitopes
        (epitope_id, protein, start_pos, end_pos, epitope_sequence,
         epitope_type, immune_response, source_database)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, insert_rows)
    db.conn.commit()

    backfill_stats = backfill_epitope_overlap(db)
    db.close()

    stats.update(backfill_stats)
    return stats


def backfill_epitope_overlap(db: DengueDB) -> dict:
    """
    Mark variant_catalog rows whose (protein, protein_position) falls inside a
    loaded epitope. Runs against variant_catalog (66k rows) rather than the
    much larger per-genome `variants` table (8M+ rows) for the same result.
    """
    cur = db.conn.cursor()

    cur.execute("SELECT protein, start_pos, end_pos, epitope_id, epitope_type FROM epitopes")
    epitopes_by_protein = defaultdict(list)
    for protein, start_pos, end_pos, epitope_id, epitope_type in cur.fetchall():
        epitopes_by_protein[protein].append((start_pos, end_pos, epitope_id, epitope_type))

    if not epitopes_by_protein:
        logger.info("No epitopes loaded; skipping variant_catalog backfill")
        return {'variant_catalog_overlaps': 0}

    cur.execute("""
        SELECT catalog_id, protein, protein_position FROM variant_catalog
        WHERE protein IS NOT NULL AND protein_position IS NOT NULL
    """)
    variant_rows = cur.fetchall()

    updates = []
    epitope_variant_counts = defaultdict(int)
    for catalog_id, protein, protein_position in variant_rows:
        hits = [(eid, etype) for start, end, eid, etype in epitopes_by_protein.get(protein, [])
                if start <= protein_position <= end]
        if hits:
            for eid, _ in hits:
                epitope_variant_counts[eid] += 1
            hit_summary = '; '.join(f"{eid}:{etype}" for eid, etype in hits)
            updates.append((hit_summary, catalog_id))

    cur.executemany(
        "UPDATE variant_catalog SET epitope_overlap = 1, epitope_hits = ? WHERE catalog_id = ?",
        updates
    )

    if epitope_variant_counts:
        cur.executemany(
            "UPDATE epitopes SET variants_in_epitope_count = ? WHERE epitope_id = ?",
            [(count, eid) for eid, count in epitope_variant_counts.items()]
        )

    db.conn.commit()
    logger.info(f"variant_catalog: flagged {len(updates)}/{len(variant_rows)} protein-coding "
                f"variants as overlapping a known epitope")
    return {'variant_catalog_overlaps': len(updates)}


def main():
    parser = argparse.ArgumentParser(description="Load IEDB epitopes and cross-reference against variant_catalog")
    parser.add_argument('--db', default='dengue_genomes.db', help='Path to SQLite database')
    parser.add_argument('--json', default=str(IEDB_JSON_PATH), help='Path to IEDB epitope JSON')
    args = parser.parse_args()

    stats = load_epitopes(args.db, Path(args.json))
    logger.info(f"Summary: {stats}")


if __name__ == "__main__":
    main()
