"""
Dengue Genome Aggregator - Fingerprint Cache Builder

Builds the per-serotype SNP fingerprint matrices used by the web app's
nearest-neighbor search (see fingerprints.py docstring for why this exists
instead of querying `variants` live).

Run once after annotate_variants.py has populated `variants`, and again any
time that table changes.

Usage:
    python build_fingerprints.py --db dengue_genomes.db
"""

import argparse
import logging
import sqlite3
import time

import references
from fingerprints import build_serotype_fingerprint

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Build SNP fingerprint cache for nearest-neighbor search")
    parser.add_argument('--db', default='dengue_genomes.db', help='Path to SQLite database')
    args = parser.parse_args()

    conn = sqlite3.connect(args.db)
    try:
        for serotype in (1, 2, 3, 4):
            ref = references.get_reference(serotype)
            t0 = time.time()
            stats = build_serotype_fingerprint(conn, serotype, ref['length'])
            logger.info(
                f"DENV-{serotype}: {stats['n_genomes']} genomes, "
                f"{stats['n_snp_variants']} SNPs fingerprinted, "
                f"{stats['n_indels_or_skipped']} indels/skipped, "
                f"{time.time() - t0:.1f}s"
            )
    finally:
        conn.close()


if __name__ == "__main__":
    main()
