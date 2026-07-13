"""
Dengue Genome Aggregator - Variant Catalog Builder

Aggregates the per-genome `variants` table (one row per variant instance per
genome) into `variant_catalog` (one row per unique variant, the browsable
"ClinVar-style" layer): occurrence count and frequency across all annotated
genomes of that serotype, date range, distinct countries, and distinct
lineages carrying it.

Run after annotate_variants.py has populated `variants`.

Usage:
    python build_variant_catalog.py --db dengue_genomes.db
"""

import argparse
import logging
import time
import uuid

from db_schema import DengueDB

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def build_variant_catalog(db_path: str) -> dict:
    db = DengueDB(db_path)
    cur = db.conn.cursor()

    # Denominator: how many genomes of each serotype actually got annotated
    # (i.e. have at least the possibility of appearing in `variants`).
    cur.execute("""
        SELECT serotype, COUNT(DISTINCT genome_id)
        FROM genomes
        WHERE genome_id IN (SELECT DISTINCT genome_id FROM variants)
        GROUP BY serotype
    """)
    genome_counts = {row[0]: row[1] for row in cur.fetchall()}
    logger.info(f"Annotated genome counts by serotype: {genome_counts}")

    logger.info("Aggregating variants (this scans the full variants table)...")
    t0 = time.time()
    cur.execute("""
        SELECT
            g.serotype,
            v.position,
            v.ref_allele,
            v.alt_allele,
            MAX(v.variant_type)     AS variant_type,
            MAX(v.protein)          AS protein,
            MAX(v.protein_position) AS protein_position,
            MAX(v.codon_change)     AS codon_change,
            MAX(v.aa_change)        AS aa_change,
            MAX(v.impact)           AS impact,
            MAX(v.is_synonymous)    AS is_synonymous,
            MAX(v.is_frameshift)    AS is_frameshift,
            MAX(v.is_stop_gained)   AS is_stop_gained,
            MAX(v.is_stop_lost)     AS is_stop_lost,
            COUNT(*)                AS occurrence_count,
            MIN(g.collection_date)  AS first_observed_date,
            MAX(g.collection_date)  AS last_observed_date,
            GROUP_CONCAT(DISTINCT g.country) AS countries,
            COUNT(DISTINCT g.country)        AS country_count,
            GROUP_CONCAT(DISTINCT g.lineage) AS lineages
        FROM variants v
        JOIN genomes g ON v.genome_id = g.genome_id
        WHERE g.serotype IS NOT NULL
        GROUP BY g.serotype, v.position, v.ref_allele, v.alt_allele
    """)
    rows = cur.fetchall()
    logger.info(f"Aggregated {len(rows)} unique variants in {time.time() - t0:.1f}s")

    cur.execute("DELETE FROM variant_catalog")

    insert_rows = []
    for r in rows:
        (serotype, position, ref_allele, alt_allele, variant_type, protein,
         protein_position, codon_change, aa_change, impact, is_synonymous,
         is_frameshift, is_stop_gained, is_stop_lost, occurrence_count,
         first_observed_date, last_observed_date, countries, country_count,
         lineages) = r

        denom = genome_counts.get(serotype, 0)
        frequency = (occurrence_count / denom) if denom else None

        insert_rows.append((
            f"VCAT_{uuid.uuid4().hex[:12]}",
            serotype, position, ref_allele, alt_allele, variant_type,
            protein, protein_position, codon_change, aa_change, impact,
            is_synonymous, is_frameshift, is_stop_gained, is_stop_lost,
            occurrence_count, denom, frequency,
            first_observed_date, last_observed_date,
            countries, country_count, lineages,
        ))

    db.conn.executemany("""
        INSERT INTO variant_catalog
        (catalog_id, serotype, position, ref_allele, alt_allele, variant_type,
         protein, protein_position, codon_change, aa_change, impact,
         is_synonymous, is_frameshift, is_stop_gained, is_stop_lost,
         occurrence_count, genome_count_serotype, frequency,
         first_observed_date, last_observed_date,
         countries, country_count, lineages)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, insert_rows)
    db.commit()

    stats = {'unique_variants': len(insert_rows), 'genome_counts_by_serotype': genome_counts}
    db.close()
    return stats


def main():
    parser = argparse.ArgumentParser(description="Build the aggregated dengue variant catalog")
    parser.add_argument('--db', default='dengue_genomes.db', help='Path to SQLite database')
    args = parser.parse_args()

    stats = build_variant_catalog(args.db)
    logger.info(f"Summary: {stats}")


if __name__ == "__main__":
    main()
