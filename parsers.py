"""
Dengue Genome Aggregator - Format Parsers
Parses FASTA, GenBank, JSON formats and normalizes to genome records
"""

import logging
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Generator
from dataclasses import dataclass
from datetime import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class GenomeRecord:
    """Standardized genome record."""
    accession: str
    serotype: Optional[int]
    lineage: Optional[str]
    genotype: Optional[str]
    clade: Optional[str]
    sequence: str
    length: int
    collection_date: Optional[str]
    collection_year: Optional[int]
    country: Optional[str]
    region: Optional[str]
    host: str = "Homo sapiens"
    source_database: Optional[str] = None
    pubmed_id: Optional[int] = None
    ncbi_taxid: int = 12637  # Dengue virus
    
    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            'accession': self.accession,
            'serotype': self.serotype,
            'lineage': self.lineage,
            'genotype': self.genotype,
            'clade': self.clade,
            'sequence': self.sequence,
            'length': self.length,
            'collection_date': self.collection_date,
            'collection_year': self.collection_year,
            'country': self.country,
            'region': self.region,
            'host': self.host,
            'source_database': self.source_database,
            'pubmed_id': self.pubmed_id,
            'ncbi_taxid': self.ncbi_taxid
        }


class FASTAParser:
    """Parse FASTA format sequences."""
    
    SEROTYPE_PATTERNS = {
        1: [r'denv1', r'denv-1', r'dengue\s*1', r'serotype\s*1', r'\bDV1\b'],
        2: [r'denv2', r'denv-2', r'dengue\s*2', r'serotype\s*2', r'\bDV2\b'],
        3: [r'denv3', r'denv-3', r'dengue\s*3', r'serotype\s*3', r'\bDV3\b'],
        4: [r'denv4', r'denv-4', r'dengue\s*4', r'serotype\s*4', r'\bDV4\b']
    }
    
    LINEAGE_PATTERNS = {
        'Asian I': [r'Asian\s*I', r'Asian-I', r'Asian1'],
        'American': [r'American', r'Americas'],
        'Cosmopolitan': [r'Cosmopolitan'],
        'Asian II': [r'Asian\s*II', r'Asian-II'],
        'African': [r'African', r'West African']
    }
    
    @classmethod
    def parse_file(cls, fasta_path: str) -> Generator[GenomeRecord, None, None]:
        """
        Parse FASTA file and yield GenomeRecord objects.
        
        Args:
            fasta_path: Path to FASTA file
        
        Yields:
            GenomeRecord objects
        """
        current_header = None
        current_seq = []
        
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('>'):
                        # Process previous record
                        if current_header and current_seq:
                            record = cls.parse_record(current_header, ''.join(current_seq))
                            if record:
                                yield record
                        
                        current_header = line[1:]  # Remove '>'
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Process last record
                if current_header and current_seq:
                    record = cls.parse_record(current_header, ''.join(current_seq))
                    if record:
                        yield record
        
        except IOError as e:
            logger.error(f"Error reading FASTA file {fasta_path}: {e}")
    
    @classmethod
    def parse_record(cls, header: str, sequence: str) -> Optional[GenomeRecord]:
        """
        Parse a single FASTA record.
        
        Args:
            header: FASTA header (after '>');
            sequence: Sequence string
        
        Returns:
            GenomeRecord or None if invalid
        """
        # Validate sequence
        if not sequence or len(sequence) < 100:
            logger.warning(f"Skipping short sequence: {header[:50]}")
            return None
        
        # Extract accession (first token, usually)
        tokens = header.split()
        accession = tokens[0] if tokens else None
        
        if not accession:
            logger.warning(f"Could not extract accession from: {header}")
            return None
        
        # Detect serotype
        serotype = cls._detect_serotype(header)
        
        # Detect lineage
        lineage = cls._detect_lineage(header)
        
        # Extract collection date
        collection_date, collection_year = cls._extract_date(header)
        
        # Extract country
        country = cls._extract_country(header)
        
        # Detect source database
        source_db = cls._detect_source_database(header)
        
        record = GenomeRecord(
            accession=accession,
            serotype=serotype,
            lineage=lineage,
            genotype=None,  # To be filled by other methods
            clade=None,
            sequence=sequence,
            length=len(sequence),
            collection_date=collection_date,
            collection_year=collection_year,
            country=country,
            region=None,  # Not easily extractable from FASTA header
            source_database=source_db
        )
        
        return record
    
    @classmethod
    def _detect_serotype(cls, header: str) -> Optional[int]:
        """Detect serotype from header text."""
        header_lower = header.lower()
        
        for serotype, patterns in cls.SEROTYPE_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, header_lower, re.IGNORECASE):
                    return serotype
        
        return None
    
    @classmethod
    def _detect_lineage(cls, header: str) -> Optional[str]:
        """Detect lineage from header text."""
        for lineage, patterns in cls.LINEAGE_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, header, re.IGNORECASE):
                    return lineage
        
        return None
    
    @classmethod
    def _extract_date(cls, header: str) -> Tuple[Optional[str], Optional[int]]:
        """Extract collection date from header."""
        # Look for YYYY-MM-DD pattern
        date_match = re.search(r'(\d{4})-(\d{2})-(\d{2})', header)
        if date_match:
            date_str = date_match.group(0)
            year = int(date_match.group(1))
            return date_str, year
        
        # Look for YYYY only
        year_match = re.search(r'\b(19\d{2}|20\d{2})\b', header)
        if year_match:
            year = int(year_match.group(1))
            return None, year
        
        return None, None
    
    @classmethod
    def _extract_country(cls, header: str) -> Optional[str]:
        """Extract country from header (naive pattern matching)."""
        # Common country names in FASTA headers
        countries = {
            'Sri Lanka': [r'sri\s*lanka', r'srilanka'],
            'Brazil': [r'brazil'],
            'Singapore': [r'singapore'],
            'Thailand': [r'thailand'],
            'Mexico': [r'mexico'],
            'USA': [r'\bUSA\b', r'\bUS\b'],
            'Indonesia': [r'indonesia'],
            'Philippines': [r'philippines'],
            'Vietnam': [r'vietnam'],
            'India': [r'india'],
            'Australia': [r'australia'],
        }
        
        for country, patterns in countries.items():
            for pattern in patterns:
                if re.search(pattern, header, re.IGNORECASE):
                    return country
        
        return None
    
    @classmethod
    def _detect_source_database(cls, header: str) -> Optional[str]:
        """Detect source database from header."""
        if re.search(r'(NC_|NZ_|AC_|NG_|NT_|NW_|XM_|XR_|YP_)', header):
            return 'NCBI'
        elif re.search(r'nextstrain', header, re.IGNORECASE):
            return 'Nextstrain'
        elif re.search(r'vipr', header, re.IGNORECASE):
            return 'ViPR'
        
        return None


class GenBankParser:
    """Parse GenBank flat file format."""
    
    @classmethod
    def parse_file(cls, genbank_path: str) -> Generator[GenomeRecord, None, None]:
        """
        Parse GenBank file and yield GenomeRecord objects.
        
        Args:
            genbank_path: Path to GenBank file
        
        Yields:
            GenomeRecord objects
        """
        try:
            from Bio import SeqIO
        except ImportError:
            logger.error("BioPython not installed. Install with: pip install biopython")
            return
        
        try:
            for record in SeqIO.parse(genbank_path, "genbank"):
                genome_record = cls.parse_record(record)
                if genome_record:
                    yield genome_record
        except Exception as e:
            logger.error(f"Error parsing GenBank file {genbank_path}: {e}")
    
    @classmethod
    def parse_record(cls, record) -> Optional[GenomeRecord]:
        """Parse a single GenBank record (BioPython SeqRecord)."""
        try:
            accession = record.id or record.name
            sequence = str(record.seq)
            
            if not sequence or len(sequence) < 100:
                return None
            
            # Extract annotations
            description = record.description or ""
            
            # Detect serotype
            serotype = cls._detect_serotype(description)
            
            # Extract from qualifiers
            collection_date = None
            collection_year = None
            country = None
            
            for feature in record.features:
                if feature.type == 'source':
                    qualifiers = feature.qualifiers
                    if 'collection_date' in qualifiers:
                        collection_date = qualifiers['collection_date'][0]
                        # Extract year
                        year_match = re.search(r'(\d{4})', collection_date)
                        if year_match:
                            collection_year = int(year_match.group(1))
                    
                    if 'country' in qualifiers:
                        country = qualifiers['country'][0]
            
            record_obj = GenomeRecord(
                accession=accession,
                serotype=serotype,
                lineage=None,
                genotype=None,
                clade=None,
                sequence=sequence,
                length=len(sequence),
                collection_date=collection_date,
                collection_year=collection_year,
                country=country,
                region=None,
                source_database='NCBI GenBank'
            )
            
            return record_obj
        
        except Exception as e:
            logger.warning(f"Error parsing GenBank record: {e}")
            return None
    
    @classmethod
    def _detect_serotype(cls, text: str) -> Optional[int]:
        """Detect serotype from description."""
        text_lower = text.lower()
        
        for i in range(1, 5):
            patterns = [
                rf'denv[- ]?{i}',
                rf'dengue\s*{i}',
                rf'serotype\s*{i}'
            ]
            for pattern in patterns:
                if re.search(pattern, text_lower):
                    return i
        
        return None


class JSONMetadataParser:
    """Parse JSON metadata (from Nextstrain, ViPR, etc.)."""
    
    @classmethod
    def parse_nextstrain_metadata(cls, metadata_dict: dict) -> Optional[GenomeRecord]:
        """Parse Nextstrain metadata JSON."""
        # Nextstrain format varies; typical structure:
        genome = GenomeRecord(
            accession=metadata_dict.get('strain') or metadata_dict.get('name'),
            serotype=None,
            lineage=None,
            genotype=metadata_dict.get('genotype'),
            clade=metadata_dict.get('clade') or metadata_dict.get('clade_membership'),
            sequence=None,  # Sequence provided separately
            length=0,
            collection_date=metadata_dict.get('date'),
            collection_year=metadata_dict.get('year'),
            country=metadata_dict.get('country'),
            region=metadata_dict.get('region'),
            source_database='Nextstrain'
        )
        
        # Try to extract year from date
        if genome.collection_date and not genome.collection_year:
            year_match = re.search(r'(\d{4})', genome.collection_date)
            if year_match:
                genome.collection_year = int(year_match.group(1))
        
        return genome
    
    @classmethod
    def parse_vipr_metadata(cls, metadata_dict: dict) -> Optional[GenomeRecord]:
        """Parse ViPR metadata JSON."""
        genome = GenomeRecord(
            accession=metadata_dict.get('accession') or metadata_dict.get('sequence_id'),
            serotype=None,
            lineage=None,
            genotype=None,
            clade=None,
            sequence=None,
            length=0,
            collection_date=metadata_dict.get('collection_date'),
            collection_year=metadata_dict.get('year'),
            country=metadata_dict.get('country'),
            region=None,
            source_database='ViPR'
        )
        
        return genome


class QualityControl:
    """Quality control checks for parsed genomes."""
    
    DENV_EXPECTED_LENGTH_RANGE = (10500, 10800)  # DENV genomes are ~10.7 kb
    MIN_SEQUENCE_IDENTITY = 0.85  # Minimum homology to dengue reference
    
    @staticmethod
    def check_sequence_length(genome: GenomeRecord) -> Tuple[bool, str]:
        """Check if sequence length is reasonable."""
        min_len, max_len = QualityControl.DENV_EXPECTED_LENGTH_RANGE
        
        if len(genome.sequence) < min_len:
            return False, f"Sequence too short: {len(genome.sequence)} bp (expected > {min_len})"
        elif len(genome.sequence) > max_len:
            return False, f"Sequence too long: {len(genome.sequence)} bp (expected < {max_len})"
        
        return True, "Length OK"
    
    @staticmethod
    def check_nucleotide_composition(genome: GenomeRecord) -> Tuple[bool, str]:
        """Check if nucleotide composition is reasonable."""
        seq = genome.sequence.upper()
        
        # Count valid nucleotides
        valid = sum(1 for nt in seq if nt in 'ATCGN')
        
        if valid / len(seq) < 0.95:
            return False, f"Too many ambiguous bases ({(1 - valid/len(seq))*100:.1f}%)"
        
        return True, "Nucleotide composition OK"
    
    @staticmethod
    def check_gc_content(genome: GenomeRecord) -> Tuple[bool, str]:
        """Check if GC content is reasonable for dengue."""
        seq = genome.sequence.upper()
        gc_count = seq.count('G') + seq.count('C')
        gc_content = gc_count / len(seq) if seq else 0
        
        # Dengue GC content typically 45-50%
        if gc_content < 0.40 or gc_content > 0.55:
            return False, f"Unusual GC content: {gc_content*100:.1f}%"
        
        return True, f"GC content OK ({gc_content*100:.1f}%)"
    
    @staticmethod
    def validate_genome(genome: GenomeRecord) -> Tuple[bool, List[str]]:
        """Run all QC checks."""
        issues = []
        
        checks = [
            QualityControl.check_sequence_length,
            QualityControl.check_nucleotide_composition,
            QualityControl.check_gc_content
        ]
        
        for check in checks:
            passed, message = check(genome)
            if not passed:
                issues.append(message)
        
        return len(issues) == 0, issues


if __name__ == "__main__":
    # Test parser
    print("Dengue Genome Parsers loaded successfully")
