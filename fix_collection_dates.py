"""
Dengue Genome Aggregator - Collection Date Cleanup (one-time migration)

`genomes.collection_date` stores partial-precision Nextstrain dates verbatim
(e.g. '1944-XX-XX', 'XXXX-XX-XX'). Because SQLite has no real DATE type,
these are compared as plain strings, and 'X' (0x58) sorts higher than any
digit -- so a query like `collection_date > date('now', '-6 months')` treats
every genome with a completely unknown date as if it were the most recent
sample, silently corrupting "recent outbreak" queries.

This migration:
  1. Adds a `collection_date_raw` column (if missing) and backfills it with
     the original, unmodified value -- no data is lost.
  2. Sets `collection_date` to NULL for any row that isn't a complete, valid
     YYYY-MM-DD date, so date-range SQL behaves correctly. `collection_year`
     (already populated independently) remains available for year-level
     filtering on partial dates.

New inserts no longer need this script -- db_schema.DengueDB.insert_genome
now does this normalization at insert time.

Usage:
    python fix_collection_dates.py --db dengue_genomes.db
"""

import argparse
import logging

from db_schema import DengueDB, normalize_collection_date

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def fix_collection_dates(db_path: str) -> dict:
    db = DengueDB(db_path)
    cur = db.conn.cursor()

    cur.execute("PRAGMA table_info(genomes)")
    columns = {row[1] for row in cur.fetchall()}
    if 'collection_date_raw' not in columns:
        logger.info("Adding collection_date_raw column...")
        cur.execute("ALTER TABLE genomes ADD COLUMN collection_date_raw VARCHAR(20)")
        db.conn.commit()

    cur.execute("""
        UPDATE genomes SET collection_date_raw = collection_date
        WHERE collection_date_raw IS NULL AND collection_date IS NOT NULL
    """)
    backfilled = cur.rowcount
    db.conn.commit()
    logger.info(f"Backfilled collection_date_raw for {backfilled} rows")

    cur.execute("SELECT genome_id, collection_date FROM genomes WHERE collection_date IS NOT NULL")
    rows = cur.fetchall()

    to_null = [(row['genome_id'],) for row in rows
               if normalize_collection_date(row['collection_date']) is None]

    logger.info(f"Nulling out {len(to_null)} of {len(rows)} non-empty collection_date values "
                f"(partial/unknown precision; raw value preserved in collection_date_raw)")
    cur.executemany("UPDATE genomes SET collection_date = NULL WHERE genome_id = ?", to_null)
    db.conn.commit()

    db.close()
    return {'checked': len(rows), 'nulled': len(to_null), 'raw_backfilled': backfilled}


def main():
    parser = argparse.ArgumentParser(description="Clean up partial-precision collection_date values")
    parser.add_argument('--db', default='dengue_genomes.db', help='Path to SQLite database')
    args = parser.parse_args()

    stats = fix_collection_dates(args.db)
    logger.info(f"Summary: {stats}")


if __name__ == "__main__":
    main()
