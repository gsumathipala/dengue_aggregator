"""
Dengue Genome Aggregator - SNP Fingerprint Cache

Precomputes, per serotype, a compact SNP profile for every annotated genome
in `variants`: one byte per reference position, 0 = matches reference,
1-4 = alt base A/C/G/T. Indels are excluded (no fixed-width positional
encoding, and a small minority of variants) -- this is a SNP-similarity
fingerprint for nearest-neighbor search, not a full variant replay.

Why this exists: a live query against the raw `variants` table (8.07M rows)
was benchmarked during development at 11.7s for a single candidate-gathering
query, because `position` alone isn't selective enough within one serotype
(nearly every genome varies at nearly every commonly-varying site). That's a
non-starter per web request. Instead, build_fingerprints.py does one full
scan up front; the web app loads the resulting small matrices and answers
nearest-neighbor queries with a single vectorized numpy comparison
(sub-100ms).
"""

import json
import logging
import sqlite3
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

CACHE_DIR = Path(__file__).parent / "fingerprint_cache"

BASE_CODE = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

_cache: dict = {}  # serotype -> {'matrix': np.ndarray, 'meta': list[dict]}


def _matrix_path(serotype: int) -> Path:
    return CACHE_DIR / f"serotype_{serotype}.npz"


def _meta_path(serotype: int) -> Path:
    return CACHE_DIR / f"serotype_{serotype}_meta.json"


def build_serotype_fingerprint(conn: sqlite3.Connection, serotype: int, ref_len: int) -> dict:
    """
    Scan `variants` + `genomes` for one serotype and write its fingerprint
    matrix + genome metadata index to fingerprint_cache/.

    Returns {'n_genomes', 'n_snp_variants', 'n_indels_or_skipped'}.
    """
    cur = conn.cursor()
    cur.execute(
        "SELECT genome_id, accession, country, collection_date, lineage, genotype "
        "FROM genomes WHERE serotype = ? ORDER BY genome_id",
        (serotype,),
    )
    genome_rows = cur.fetchall()
    genome_ids = [r[0] for r in genome_rows]
    index = {gid: i for i, gid in enumerate(genome_ids)}

    n = len(genome_ids)
    matrix = np.zeros((n, ref_len), dtype=np.uint8)

    cur.execute(
        """
        SELECT v.genome_id, v.position, v.alt_allele, v.variant_type
        FROM variants v
        JOIN genomes g ON v.genome_id = g.genome_id
        WHERE g.serotype = ?
        """,
        (serotype,),
    )

    n_snp = 0
    n_skipped = 0
    for genome_id, position, alt_allele, variant_type in cur:
        row = index.get(genome_id)
        if row is None:
            n_skipped += 1
            continue
        code = BASE_CODE.get(alt_allele) if variant_type == 'SNP' else None
        if code is None or not (1 <= position <= ref_len):
            n_skipped += 1
            continue
        matrix[row, position - 1] = code
        n_snp += 1

    CACHE_DIR.mkdir(exist_ok=True)
    np.savez_compressed(_matrix_path(serotype), matrix=matrix)

    meta = [
        {
            'genome_id': gid,
            'accession': acc,
            'country': country,
            'collection_date': date,
            'lineage': lineage,
            'genotype': genotype,
            'variant_count': int(np.count_nonzero(matrix[i])),
        }
        for i, (gid, acc, country, date, lineage, genotype) in enumerate(genome_rows)
    ]
    _meta_path(serotype).write_text(json.dumps(meta), encoding='utf-8')

    return {'n_genomes': n, 'n_snp_variants': n_snp, 'n_indels_or_skipped': n_skipped}


def cache_exists(serotype: int) -> bool:
    return _matrix_path(serotype).exists() and _meta_path(serotype).exists()


def all_caches_exist() -> bool:
    return all(cache_exists(s) for s in (1, 2, 3, 4))


def load_fingerprints(serotype: int) -> dict:
    """Load (and memoize in-process) the cached matrix + metadata for one serotype."""
    if serotype in _cache:
        return _cache[serotype]

    if not cache_exists(serotype):
        raise FileNotFoundError(
            f"No fingerprint cache for serotype {serotype}. Run build_fingerprints.py first."
        )

    with np.load(_matrix_path(serotype)) as data:
        matrix = data['matrix']
    meta = json.loads(_meta_path(serotype).read_text(encoding='utf-8'))

    result = {'matrix': matrix, 'meta': meta}
    _cache[serotype] = result
    return result


def nearest_neighbors(serotype: int, query_vec: np.ndarray, covered_mask: np.ndarray, top_n: int = 10) -> list:
    """
    Find the top_n genomes (of `serotype`) whose SNP profile best explains the
    query's own called variants.

    Args:
        query_vec: uint8 array, same length as the reference, coded like the
                   cached matrix (0 = ref, 1-4 = alt base), built from the
                   query's own called SNPs.
        covered_mask: boolean array, same length, True where the query
                      actually has aligned sequence data (its coverage span).
        top_n: how many candidates to return.

    Returns:
        List of dicts (richest match first): genome metadata plus
        shared_variants, query_variant_count, candidate_variant_count, and
        similarity_pct = 100 * shared / query_variant_count. Similarity is
        scored only over positions where the *query* has a variant, so it
        answers "of the mutations your sequence carries, what fraction does
        this database genome also carry" -- meaningful for partial sequences
        regardless of how much of the genome they cover.
    """
    data = load_fingerprints(serotype)
    matrix = data['matrix']
    meta = data['meta']

    query_variant_mask = covered_mask & (query_vec != 0)
    query_variant_count = int(np.count_nonzero(query_variant_mask))

    if query_variant_count == 0 or matrix.shape[0] == 0:
        return []

    sub_matrix = matrix[:, query_variant_mask]
    sub_query = query_vec[query_variant_mask]
    shared = (sub_matrix == sub_query).sum(axis=1)

    order = np.argsort(-shared)[:top_n]

    results = []
    for i in order:
        shared_count = int(shared[i])
        if shared_count == 0:
            continue
        m = meta[i]
        results.append({
            'accession': m['accession'],
            'country': m['country'],
            'collection_date': m['collection_date'],
            'lineage': m['lineage'],
            'genotype': m['genotype'],
            'shared_variants': shared_count,
            'query_variant_count': query_variant_count,
            'candidate_variant_count': m['variant_count'],
            'similarity_pct': round(100 * shared_count / query_variant_count, 1),
        })
    return results
