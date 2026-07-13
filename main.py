"""
Dengue Virus Genome Aggregator - Main Entry Point
Complete pipeline orchestration with CLI interface
"""

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional

from db_schema import DengueDB
from downloaders import (
    NCBIDatasetsDownloader, ViPRDownloader, 
    NextstrainDownloader, download_all_sources
)
from parsers import FASTAParser, QualityControl
from variant_annotator import VariantAnnotator, SequenceAligner

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DengueAggregatorPipeline:
    """Main dengue aggregator pipeline."""
    
    def __init__(self, db_path: str = "dengue_genomes.db", data_dir: str = "dengue_data"):
        self.db_path = db_path
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        self.db = DengueDB(db_path)
        self.stats = {
            'downloaded': 0,
            'parsed': 0,
            'passed_qc': 0,
            'failed_qc': 0,
            'inserted': 0,
            'variants_called': 0,
            'errors': []
        }
    
    def run_full_pipeline(self, download: bool = True, parse: bool = True, 
                         annotate: bool = True) -> dict:
        """Run complete pipeline."""
        logger.info("=" * 70)
        logger.info("DENGUE VIRUS GENOME AGGREGATOR")
        logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 70)
        
        if download:
            self._download_all_sources()
        
        if parse:
            self._parse_all_files()
        
        if annotate:
            self._annotate_variants()
        
        self._print_summary()
        return self.stats
    
    def _download_all_sources(self):
        """Download from all sources."""
        logger.info("\n[PHASE 1] Downloading from all data sources...")
        stats = download_all_sources(str(self.data_dir))
        self.stats['downloaded'] = sum(stats.values())
    
    def _parse_all_files(self):
        """Parse downloaded files and insert into database."""
        logger.info("\n[PHASE 2] Parsing downloaded files...")

        fasta_files = list(self.data_dir.glob("*.fasta"))

        for fasta_file in fasta_files:
            logger.info(f"Parsing {fasta_file.name}...")

            # Extract serotype from filename
            serotype = self._extract_serotype(fasta_file.name)

            # Nextstrain FASTA headers are bare accessions with no embedded
            # metadata, so pull country/date/lineage from the sibling TSV
            # (saved alongside the FASTA as "<name>_metadata.tsv") instead.
            metadata_path = self.data_dir / f"{fasta_file.stem}_metadata.tsv"
            metadata = self._load_metadata_tsv(metadata_path) if metadata_path.exists() else {}
            if metadata:
                logger.info(f"  Using metadata from {metadata_path.name} ({len(metadata)} records)")

            source_database = self._detect_source_from_filename(fasta_file.name)

            record_count = 0
            for genome in FASTAParser.parse_file(str(fasta_file)):
                if serotype:
                    genome.serotype = serotype
                if source_database:
                    genome.source_database = source_database

                self._apply_metadata(genome, metadata.get(genome.accession))

                # QC check
                passed_qc, issues = QualityControl.validate_genome(genome)

                if passed_qc:
                    try:
                        genome_id = self.db.insert_genome(genome.to_dict())
                        self.stats['passed_qc'] += 1
                        self.stats['inserted'] += 1
                        record_count += 1
                    except Exception as e:
                        logger.error(f"Error inserting {genome.accession}: {e}")
                        self.stats['errors'].append(str(e))
                else:
                    self.stats['failed_qc'] += 1
                    if record_count < 5:  # Log first 5 QC failures
                        logger.debug(f"QC failed for {genome.accession}: {issues}")

            logger.info(f"  [OK] Inserted {record_count} genomes from {fasta_file.name}")
            self.stats['parsed'] += record_count

    def _extract_serotype(self, filename: str) -> int:
        """Extract serotype from filename."""
        for i in range(1, 5):
            if f"denv{i}" in filename.lower():
                return i
        return None

    @staticmethod
    def _detect_source_from_filename(filename: str) -> Optional[str]:
        """Map a downloaded FASTA filename back to the source that produced it."""
        name = filename.lower()
        if name.startswith('nextstrain'):
            return 'Nextstrain'
        if name.startswith('bvbrc'):
            return 'BV-BRC'
        if name.startswith('ncbi'):
            return 'NCBI'
        return None

    @staticmethod
    def _load_metadata_tsv(path: Path) -> dict:
        """Load a Nextstrain-style metadata TSV into a dict keyed by accession."""
        metadata = {}
        try:
            with open(path, newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    accession = row.get('accession')
                    if accession:
                        metadata[accession] = row
        except (IOError, csv.Error) as e:
            logger.warning(f"Could not read metadata file {path}: {e}")
        return metadata

    @staticmethod
    def _apply_metadata(genome, meta_row: Optional[dict]):
        """Overlay metadata TSV fields onto a parsed GenomeRecord, when available."""
        if not meta_row:
            return

        genome.country = meta_row.get('country') or genome.country
        genome.collection_date = meta_row.get('date') or genome.collection_date
        genome.genotype = meta_row.get('genotype') or genome.genotype
        genome.lineage = meta_row.get('major_lineage') or genome.genotype or genome.lineage
        genome.clade = meta_row.get('minor_lineage') or genome.clade
        if meta_row.get('region'):
            genome.region = meta_row['region']

        if genome.collection_date and not genome.collection_year:
            year_match = re.search(r'(\d{4})', genome.collection_date)
            if year_match:
                genome.collection_year = int(year_match.group(1))
    
    def _annotate_variants(self):
        """Annotate variants (placeholder - requires reference sequences)."""
        logger.info("\n[PHASE 3] Annotating variants...")
        logger.info("  Note: Variant annotation requires reference genomes.")
        logger.info("  Download from NCBI RefSeq: NC_001474 (DENV-1), NC_001438 (DENV-2), etc.")
        logger.info("  Skipping for now - can be integrated after reference setup")
    
    def _print_summary(self):
        """Print execution summary."""
        logger.info("\n" + "=" * 70)
        logger.info("EXECUTION SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Downloaded:        {self.stats['downloaded']:>6} sequences")
        logger.info(f"Parsed:            {self.stats['parsed']:>6} sequences")
        logger.info(f"Passed QC:         {self.stats['passed_qc']:>6} sequences")
        logger.info(f"Failed QC:         {self.stats['failed_qc']:>6} sequences")
        logger.info(f"Inserted to DB:    {self.stats['inserted']:>6} sequences")
        logger.info(f"Errors:            {len(self.stats['errors']):>6} issues")
        logger.info("=" * 70)
        logger.info(f"Database: {self.db_path}")
        logger.info(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 70 + "\n")
    
    def query_genomes(self, serotype: int = None, country: str = None):
        """Query genomes from database."""
        filters = {}
        if serotype:
            filters['serotype'] = serotype
        if country:
            filters['country'] = country
        
        results = self.db.query_genomes(filters)
        logger.info(f"Found {len(results)} matching genomes")
        
        for row in results[:10]:  # Show first 10
            logger.info(f"  {row['accession']}: DENV-{row['serotype']} from {row['country']} ({row['collection_date']})")
        
        return results
    
    def close(self):
        """Close database connection."""
        self.db.close()


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Dengue Virus Genome Aggregator - Download, parse, and annotate all DENV sequences"
    )
    
    parser.add_argument('--db', default='dengue_genomes.db',
                       help='Path to SQLite database (default: dengue_genomes.db)')
    parser.add_argument('--data-dir', default='dengue_data',
                       help='Directory for downloaded data (default: dengue_data)')
    parser.add_argument('--skip-download', action='store_true',
                       help='Skip download phase (use existing files)')
    parser.add_argument('--skip-parse', action='store_true',
                       help='Skip parse phase')
    parser.add_argument('--skip-annotate', action='store_true',
                       help='Skip variant annotation phase')
    parser.add_argument('--query-serotype', type=int,
                       help='Query genomes by serotype (1-4)')
    parser.add_argument('--query-country', type=str,
                       help='Query genomes by country')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize pipeline
    pipeline = DengueAggregatorPipeline(db_path=args.db, data_dir=args.data_dir)
    
    try:
        # Run pipeline
        pipeline.run_full_pipeline(
            download=not args.skip_download,
            parse=not args.skip_parse,
            annotate=not args.skip_annotate
        )
        
        # Optional queries
        if args.query_serotype or args.query_country:
            logger.info("\nRunning queries...")
            pipeline.query_genomes(
                serotype=args.query_serotype,
                country=args.query_country
            )
    
    finally:
        pipeline.close()


if __name__ == "__main__":
    main()
