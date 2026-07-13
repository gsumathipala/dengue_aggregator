"""
Dengue Genome Aggregator - Bulk Variant Calling
Aligns every genome in the database to its serotype's reference genome and
populates the `variants` table. Run after main.py has populated `genomes`.

Usage:
    python annotate_variants.py --db dengue_genomes.db
"""

import argparse
import logging
import time

from db_schema import DengueDB
from variant_annotator import VariantAnnotator

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def annotate_all_variants(db_path: str, commit_every: int = 200, only_missing: bool = True) -> dict:
    """
    Call variants for every genome in the database and insert them.

    Args:
        db_path: Path to the SQLite database (must already have a populated
                 `genomes` table).
        commit_every: Commit after this many genomes' worth of variants.
        only_missing: Skip genomes that already have variant rows (lets this
                      be re-run safely after new genomes are added).

    Returns:
        Stats dict: genomes_processed, variants_inserted, errors.
    """
    db = DengueDB(db_path)
    stats = {'genomes_processed': 0, 'variants_inserted': 0, 'errors': 0}

    for serotype in (1, 2, 3, 4):
        cur = db.conn.cursor()
        query = "SELECT genome_id, accession, sequence FROM genomes WHERE serotype = ?"
        if only_missing:
            query += " AND genome_id NOT IN (SELECT DISTINCT genome_id FROM variants)"
        cur.execute(query, (serotype,))
        rows = cur.fetchall()

        if not rows:
            logger.info(f"DENV-{serotype}: nothing to do")
            continue

        annotator = VariantAnnotator(serotype=serotype)
        logger.info(f"DENV-{serotype}: annotating {len(rows)} genomes...")

        t_start = time.time()
        since_commit = 0

        for i, row in enumerate(rows):
            genome_id, accession, sequence = row['genome_id'], row['accession'], row['sequence']
            if not sequence:
                continue

            try:
                variants = annotator.call_variants(sequence)
                db.insert_variants_bulk(genome_id, [v.to_dict() for v in variants])
                stats['variants_inserted'] += len(variants)
            except Exception as e:
                logger.error(f"Error annotating {accession}: {e}")
                stats['errors'] += 1
                continue

            stats['genomes_processed'] += 1
            since_commit += 1
            if since_commit >= commit_every:
                db.commit()
                since_commit = 0

            if (i + 1) % 1000 == 0:
                elapsed = time.time() - t_start
                rate = (i + 1) / elapsed
                logger.info(f"  DENV-{serotype}: {i + 1}/{len(rows)} genomes ({rate:.1f}/s)")

        db.commit()
        elapsed = time.time() - t_start
        logger.info(f"DENV-{serotype}: done in {elapsed:.1f}s")

    db.close()
    return stats


def main():
    parser = argparse.ArgumentParser(description="Bulk DENV variant calling against serotype references")
    parser.add_argument('--db', default='dengue_genomes.db', help='Path to SQLite database')
    parser.add_argument('--commit-every', type=int, default=200, help='Commit interval (genomes)')
    parser.add_argument('--all', dest='only_missing', action='store_false',
                         help='Re-annotate genomes that already have variant rows')
    args = parser.parse_args()

    stats = annotate_all_variants(args.db, commit_every=args.commit_every, only_missing=args.only_missing)
    logger.info(f"Summary: {stats}")


if __name__ == "__main__":
    main()
