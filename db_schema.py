"""
Dengue Virus Genome Aggregator - Database Schema
Initializes SQLite/PostgreSQL with all required tables and indices
"""

import sqlite3
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def normalize_collection_date(raw_date: Optional[str]) -> Optional[str]:
    """
    Return raw_date unchanged if it's a complete, valid YYYY-MM-DD date, else None.

    Sources like Nextstrain report partial-precision dates with 'X' placeholders
    for an unknown month/day (e.g. '1944-XX-XX', 'XXXX-XX-XX'). Storing those
    verbatim in `collection_date` corrupts SQL date comparisons: 'X' sorts
    higher than any digit, so 'XXXX-XX-XX' > date('now', '-6 months') is TRUE,
    making fully-unknown dates look like the most recent ones. The raw string
    is preserved separately in `collection_date_raw`; `collection_year` (parsed
    independently) still supports year-level filtering for partial dates.
    """
    if not raw_date:
        return None
    try:
        datetime.strptime(raw_date, '%Y-%m-%d')
        return raw_date
    except ValueError:
        return None


class DengueDB:
    """Database initialization and management for dengue genome aggregator."""
    
    def __init__(self, db_path: str = "dengue_genomes.db", db_type: str = "sqlite"):
        """
        Initialize database connection.
        
        Args:
            db_path: Path to SQLite database file
            db_type: 'sqlite' or 'postgresql'
        """
        self.db_path = db_path
        self.db_type = db_type
        self.conn = None
        self.cursor = None
        
        if db_type == "sqlite":
            self._init_sqlite()
        elif db_type == "postgresql":
            raise NotImplementedError("PostgreSQL support coming soon")
        else:
            raise ValueError(f"Unsupported db_type: {db_type}")
    
    def _init_sqlite(self):
        """Initialize SQLite connection and create schema."""
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.cursor = self.conn.cursor()
        logger.info(f"Connected to SQLite database: {self.db_path}")
        self._create_schema()
    
    def _create_schema(self):
        """Create all required tables and indices."""
        schema_sql = """
        -- Main genomes table
        CREATE TABLE IF NOT EXISTS genomes (
            genome_id TEXT PRIMARY KEY,
            accession VARCHAR(20) UNIQUE NOT NULL,
            serotype INTEGER CHECK(serotype IN (1, 2, 3, 4)),
            lineage VARCHAR(50),
            genotype VARCHAR(100),
            clade VARCHAR(50),
            length INTEGER,
            sequence TEXT,
            collection_date DATE,
            collection_date_raw VARCHAR(20),
            collection_year INTEGER,
            country VARCHAR(100),
            region VARCHAR(100),
            host TEXT DEFAULT 'Homo sapiens',
            source_database VARCHAR(50),
            pubmed_id INTEGER,
            ncbi_taxid INTEGER,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        
        CREATE INDEX IF NOT EXISTS idx_genomes_serotype ON genomes(serotype);
        CREATE INDEX IF NOT EXISTS idx_genomes_country ON genomes(country);
        CREATE INDEX IF NOT EXISTS idx_genomes_collection_date ON genomes(collection_date);
        CREATE INDEX IF NOT EXISTS idx_genomes_lineage ON genomes(lineage);
        CREATE INDEX IF NOT EXISTS idx_genomes_accession ON genomes(accession);
        
        -- Variants table (SNPs, indels, MNPs)
        CREATE TABLE IF NOT EXISTS variants (
            variant_id TEXT PRIMARY KEY,
            genome_id TEXT NOT NULL,
            position INTEGER NOT NULL,
            ref_allele VARCHAR(100),
            alt_allele VARCHAR(100),
            variant_type VARCHAR(20),
            annotation TEXT,
            impact VARCHAR(20),
            protein VARCHAR(50),
            protein_position INTEGER,
            codon_change VARCHAR(20),
            aa_change VARCHAR(20),
            is_synonymous BOOLEAN,
            is_frameshift BOOLEAN,
            is_stop_gained BOOLEAN,
            is_stop_lost BOOLEAN,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
        );
        
        CREATE INDEX IF NOT EXISTS idx_variants_genome ON variants(genome_id);
        CREATE INDEX IF NOT EXISTS idx_variants_position ON variants(position);
        CREATE INDEX IF NOT EXISTS idx_variants_protein ON variants(protein);
        CREATE INDEX IF NOT EXISTS idx_variants_impact ON variants(impact);
        CREATE INDEX IF NOT EXISTS idx_variants_aa_change ON variants(aa_change);
        
        -- Proteins table
        CREATE TABLE IF NOT EXISTS proteins (
            protein_id TEXT PRIMARY KEY,
            genome_id TEXT NOT NULL,
            protein_name VARCHAR(50),
            protein_description TEXT,
            start_pos INTEGER,
            end_pos INTEGER,
            length INTEGER,
            sequence TEXT,
            function TEXT,
            pdb_id VARCHAR(10),
            esm_embedding BLOB,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
        );
        
        CREATE INDEX IF NOT EXISTS idx_proteins_genome ON proteins(genome_id);
        CREATE INDEX IF NOT EXISTS idx_proteins_name ON proteins(protein_name);
        
        -- DENV proteins (standard annotation)
        CREATE TABLE IF NOT EXISTS protein_annotations (
            annotation_id TEXT PRIMARY KEY,
            protein_name VARCHAR(50),
            serotype INTEGER,
            start_pos INTEGER,
            end_pos INTEGER,
            description TEXT,
            function TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        
        -- Epitopes table (B-cell, T-cell)
        CREATE TABLE IF NOT EXISTS epitopes (
            epitope_id TEXT PRIMARY KEY,
            protein VARCHAR(50),
            start_pos INTEGER,
            end_pos INTEGER,
            epitope_sequence VARCHAR(100),
            epitope_type VARCHAR(50),
            immune_response VARCHAR(100),
            source_database VARCHAR(50),
            variants_in_epitope_count INTEGER DEFAULT 0,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        
        CREATE INDEX IF NOT EXISTS idx_epitopes_protein ON epitopes(protein);
        CREATE INDEX IF NOT EXISTS idx_epitopes_position ON epitopes(start_pos, end_pos);
        
        -- Lineage and clade assignments
        CREATE TABLE IF NOT EXISTS lineage_assignments (
            assignment_id TEXT PRIMARY KEY,
            genome_id TEXT NOT NULL,
            serotype INTEGER,
            genotype VARCHAR(50),
            clade VARCHAR(50),
            assignment_method VARCHAR(50),
            confidence REAL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
        );
        
        CREATE INDEX IF NOT EXISTS idx_lineage_genome ON lineage_assignments(genome_id);
        CREATE INDEX IF NOT EXISTS idx_lineage_genotype ON lineage_assignments(genotype);
        
        -- Metadata (flexible key-value storage)
        CREATE TABLE IF NOT EXISTS metadata (
            metadata_id TEXT PRIMARY KEY,
            genome_id TEXT NOT NULL,
            key VARCHAR(100),
            value TEXT,
            source_database VARCHAR(50),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
        );
        
        CREATE INDEX IF NOT EXISTS idx_metadata_genome ON metadata(genome_id);
        CREATE INDEX IF NOT EXISTS idx_metadata_key ON metadata(key);
        
        -- Source tracking (which databases contributed which genomes)
        CREATE TABLE IF NOT EXISTS source_tracking (
            tracking_id TEXT PRIMARY KEY,
            genome_id TEXT NOT NULL,
            source_name VARCHAR(50),
            source_accession VARCHAR(100),
            download_date DATE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
        );
        
        -- Integrity checks and data quality
        CREATE TABLE IF NOT EXISTS quality_metrics (
            metric_id TEXT PRIMARY KEY,
            genome_id TEXT NOT NULL,
            sequence_length INTEGER,
            gc_content REAL,
            n_count INTEGER,
            ambiguous_bases INTEGER,
            quality_score REAL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY(genome_id) REFERENCES genomes(genome_id)
        );
        
        CREATE INDEX IF NOT EXISTS idx_quality_genome ON quality_metrics(genome_id);
        
        -- Log table for tracking updates
        CREATE TABLE IF NOT EXISTS update_log (
            log_id INTEGER PRIMARY KEY AUTOINCREMENT,
            event_type VARCHAR(50),
            source_database VARCHAR(50),
            genomes_added INTEGER,
            genomes_updated INTEGER,
            variants_added INTEGER,
            timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );

        -- Aggregated, one-row-per-unique-variant browsable catalog (the
        -- "ClinVar-style" layer): population frequency across all annotated
        -- genomes, plus optional epitope-overlap and literature annotations.
        CREATE TABLE IF NOT EXISTS variant_catalog (
            catalog_id TEXT PRIMARY KEY,
            serotype INTEGER NOT NULL,
            position INTEGER NOT NULL,
            ref_allele VARCHAR(100),
            alt_allele VARCHAR(100),
            variant_type VARCHAR(20),
            protein VARCHAR(50),
            protein_position INTEGER,
            codon_change VARCHAR(20),
            aa_change VARCHAR(20),
            impact VARCHAR(20),
            is_synonymous BOOLEAN,
            is_frameshift BOOLEAN,
            is_stop_gained BOOLEAN,
            is_stop_lost BOOLEAN,
            occurrence_count INTEGER,
            genome_count_serotype INTEGER,
            frequency REAL,
            first_observed_date DATE,
            last_observed_date DATE,
            countries TEXT,
            country_count INTEGER,
            lineages TEXT,
            epitope_overlap BOOLEAN DEFAULT 0,
            epitope_hits TEXT,
            literature_note TEXT,
            literature_refs TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );

        CREATE UNIQUE INDEX IF NOT EXISTS idx_variant_catalog_unique
            ON variant_catalog(serotype, position, ref_allele, alt_allele);
        CREATE INDEX IF NOT EXISTS idx_variant_catalog_protein ON variant_catalog(protein, aa_change);
        CREATE INDEX IF NOT EXISTS idx_variant_catalog_impact ON variant_catalog(impact);
        CREATE INDEX IF NOT EXISTS idx_variant_catalog_frequency ON variant_catalog(frequency);
        CREATE INDEX IF NOT EXISTS idx_variant_catalog_epitope ON variant_catalog(epitope_overlap);
        """
        
        # Split and execute each statement
        for statement in schema_sql.split(';'):
            if statement.strip():
                try:
                    self.cursor.execute(statement)
                except sqlite3.OperationalError as e:
                    # Table might already exist
                    logger.debug(f"Schema statement skipped: {e}")
        
        self.conn.commit()
        logger.info("Database schema initialized successfully")
    
    def insert_genome(self, genome_data: dict) -> str:
        """
        Insert a genome record.
        
        Args:
            genome_data: Dictionary with keys: accession, serotype, lineage, 
                        length, sequence, collection_date, country, etc.
        
        Returns:
            genome_id (generated)
        """
        import uuid
        genome_id = f"DENV_{uuid.uuid4().hex[:12]}"
        
        sql = """
        INSERT INTO genomes
        (genome_id, accession, serotype, lineage, genotype, clade, length, sequence,
         collection_date, collection_date_raw, collection_year, country, region, host,
         source_database, pubmed_id, ncbi_taxid)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        raw_date = genome_data.get('collection_date')
        values = (
            genome_id,
            genome_data.get('accession'),
            genome_data.get('serotype'),
            genome_data.get('lineage'),
            genome_data.get('genotype'),
            genome_data.get('clade'),
            genome_data.get('length'),
            genome_data.get('sequence'),
            normalize_collection_date(raw_date),
            raw_date,
            genome_data.get('collection_year'),
            genome_data.get('country'),
            genome_data.get('region'),
            genome_data.get('host', 'Homo sapiens'),
            genome_data.get('source_database'),
            genome_data.get('pubmed_id'),
            genome_data.get('ncbi_taxid')
        )
        
        self.cursor.execute(sql, values)
        self.conn.commit()
        logger.info(f"Inserted genome: {genome_data.get('accession')} -> {genome_id}")
        return genome_id
    
    def insert_variant(self, variant_data: dict) -> str:
        """Insert a single variant record (commits immediately)."""
        import uuid
        variant_id = f"VAR_{uuid.uuid4().hex[:12]}"

        sql = """
        INSERT INTO variants
        (variant_id, genome_id, position, ref_allele, alt_allele,
         variant_type, impact, protein, protein_position, codon_change,
         aa_change, is_synonymous, is_frameshift, is_stop_gained, is_stop_lost)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        values = (
            variant_id,
            variant_data.get('genome_id'),
            variant_data.get('position'),
            variant_data.get('ref_allele'),
            variant_data.get('alt_allele'),
            variant_data.get('variant_type'),
            variant_data.get('impact'),
            variant_data.get('protein'),
            variant_data.get('protein_position'),
            variant_data.get('codon_change'),
            variant_data.get('aa_change'),
            variant_data.get('is_synonymous'),
            variant_data.get('is_frameshift'),
            variant_data.get('is_stop_gained'),
            variant_data.get('is_stop_lost')
        )

        self.cursor.execute(sql, values)
        self.conn.commit()
        return variant_id

    def insert_variants_bulk(self, genome_id: str, variants: list) -> int:
        """
        Insert many variant records for one genome in a single transaction
        (no per-row commit). Intended for bulk variant-calling runs where
        per-row commits would be prohibitively slow.

        Args:
            genome_id: The genome these variants belong to.
            variants: List of variant dicts (same shape as insert_variant's
                      variant_data, minus genome_id).

        Returns:
            Number of rows inserted.
        """
        import uuid

        sql = """
        INSERT INTO variants
        (variant_id, genome_id, position, ref_allele, alt_allele,
         variant_type, impact, protein, protein_position, codon_change,
         aa_change, is_synonymous, is_frameshift, is_stop_gained, is_stop_lost)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        rows = [
            (
                f"VAR_{uuid.uuid4().hex[:12]}",
                genome_id,
                v.get('position'),
                v.get('ref_allele'),
                v.get('alt_allele'),
                v.get('variant_type'),
                v.get('impact'),
                v.get('protein'),
                v.get('protein_position'),
                v.get('codon_change'),
                v.get('aa_change'),
                v.get('is_synonymous'),
                v.get('is_frameshift'),
                v.get('is_stop_gained'),
                v.get('is_stop_lost'),
            )
            for v in variants
        ]

        self.cursor.executemany(sql, rows)
        return len(rows)

    def query_genomes(self, filters: dict) -> list:
        """
        Query genomes with flexible filters.
        
        Args:
            filters: Dict with optional keys: serotype, country, lineage, 
                    after_date, before_date
        
        Returns:
            List of genome records
        """
        sql = "SELECT * FROM genomes WHERE 1=1"
        params = []
        
        if 'serotype' in filters:
            sql += " AND serotype = ?"
            params.append(filters['serotype'])
        
        if 'country' in filters:
            sql += " AND country = ?"
            params.append(filters['country'])
        
        if 'lineage' in filters:
            sql += " AND lineage = ?"
            params.append(filters['lineage'])
        
        if 'after_date' in filters:
            sql += " AND collection_date >= ?"
            params.append(filters['after_date'])
        
        if 'before_date' in filters:
            sql += " AND collection_date <= ?"
            params.append(filters['before_date'])
        
        sql += " ORDER BY collection_date DESC"

        self.cursor.execute(sql, params)
        return self.cursor.fetchall()

    def commit(self):
        """Commit the current transaction (used after batched bulk inserts)."""
        self.conn.commit()

    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()
            logger.info("Database connection closed")
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


if __name__ == "__main__":
    # Test database initialization
    db = DengueDB("dengue_test.db")
    logger.info("Database initialized successfully")
    db.close()
