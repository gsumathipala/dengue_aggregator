"""
Dengue Genome Aggregator - Multi-Source Downloaders
Fetches sequences from NCBI, BV-BRC (successor to ViPR), and Nextstrain.

NOTE (repair, 2026-07-12): The previous version of this module targeted
endpoints that no longer exist or never existed in the shape assumed here:
  - "NCBI Datasets v1" virus/taxon_tree, virus/by_accession, virus/download
    endpoints returned 0 results (wrong API shape).
  - ViPR (viprbrc.org) has been decommissioned; the site now redirects to
    BV-BRC (bv-brc.org), which returned an HTML homepage instead of JSON.
  - Nextstrain's old "charon/v2/dengue/.../sequences.fasta" URLs 404.
All three were verified live and replaced below with endpoints confirmed
to return real data:
  - NCBI: classic E-utilities (esearch/efetch against db=nuccore).
  - BV-BRC: its documented Solr/RQL Data API (genome + genome_sequence).
  - Nextstrain: the curated ingest files published to data.nextstrain.org
    (zstd-compressed FASTA/TSV), as used by nextstrain/dengue itself.
"""

import requests
import logging
import time
from pathlib import Path
from typing import List, Dict, Optional, Generator, Tuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class NCBIDatasetsDownloader:
    """Download DENV sequences from NCBI via E-utilities (esearch/efetch)."""

    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    DENV_TAXON_ID = "12637"  # Dengue virus taxon ID
    SEARCH_TERM = 'txid12637[Organism:exp] AND "complete genome"[Title]'

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize NCBI downloader.

        Args:
            api_key: Optional NCBI API key for higher rate limits
                     (10 req/s vs. 3 req/s without one).
        """
        self.api_key = api_key
        self.session = requests.Session()

    def _esearch(self, term: str, retmax: int) -> dict:
        """Run an esearch query with the history server enabled."""
        params = {
            "db": "nuccore",
            "term": term,
            "retmode": "json",
            "retmax": retmax,
            "usehistory": "y",
        }
        if self.api_key:
            params["api_key"] = self.api_key

        response = self.session.get(f"{self.EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
        response.raise_for_status()
        return response.json()["esearchresult"]

    def download_dengue_genome_list(self, max_genomes: int = 500) -> List[Dict]:
        """
        Retrieve accessions for complete DENV genomes from NCBI.

        Args:
            max_genomes: Maximum number of UIDs to retrieve.

        Returns:
            List of dicts with 'uid' and 'source_database'.
        """
        logger.info("Searching NCBI nuccore for DENV complete genomes...")
        try:
            result = self._esearch(self.SEARCH_TERM, retmax=max_genomes)
            count = int(result.get("count", 0))
            ids = result.get("idlist", [])
            logger.info(f"NCBI reports {count} matching records; retrieved {len(ids)} UIDs")
            return [{"uid": uid, "source_database": "NCBI"} for uid in ids]
        except (requests.RequestException, ValueError, KeyError) as e:
            logger.error(f"Error searching NCBI: {e}")
            return []

    def download_genome_sequence(self, accession: str) -> Optional[str]:
        """
        Download a single FASTA sequence by accession via efetch.

        Args:
            accession: GenBank accession (e.g., NC_001474)

        Returns:
            FASTA text (including header) or None if not found.
        """
        params = {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"}
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            response = self.session.get(f"{self.EUTILS_BASE}/efetch.fcgi", params=params, timeout=30)
            response.raise_for_status()
            return response.text if response.text.startswith(">") else None
        except requests.RequestException as e:
            logger.error(f"Error downloading sequence {accession}: {e}")
            return None

    def bulk_download_fasta(self, output_dir: str = "ncbi_dengue",
                             max_genomes: int = 500) -> Generator[Tuple[str, str], None, None]:
        """
        Bulk download DENV complete genomes as FASTA using the NCBI history
        server (esearch + batched efetch), respecting NCBI rate limits.

        Args:
            output_dir: Present for API-compatibility with the old signature;
                        this method only yields records, it does not write files.
            max_genomes: Cap on how many genomes to fetch (NCBI currently has
                         several thousand; capped by default to keep runtime sane).

        Yields:
            (header, sequence) tuples
        """
        try:
            search = self._esearch(self.SEARCH_TERM, retmax=max_genomes)
        except (requests.RequestException, ValueError, KeyError) as e:
            logger.error(f"Error searching NCBI: {e}")
            return

        webenv = search.get("webenv")
        query_key = search.get("querykey")
        total = min(int(search.get("count", 0)), max_genomes)

        if total > int(search.get("count", 0)):
            total = int(search.get("count", 0))
        if int(search.get("count", 0)) > max_genomes:
            logger.info(f"NCBI has {search.get('count')} matching genomes; "
                        f"downloading the first {max_genomes} (increase max_genomes for more)")

        logger.info(f"Fetching {total} DENV genomes from NCBI via E-utilities...")

        batch_size = 200
        for start in range(0, total, batch_size):
            params = {
                "db": "nuccore",
                "rettype": "fasta",
                "retmode": "text",
                "retstart": start,
                "retmax": min(batch_size, total - start),
                "webenv": webenv,
                "query_key": query_key,
            }
            if self.api_key:
                params["api_key"] = self.api_key

            try:
                response = self.session.get(f"{self.EUTILS_BASE}/efetch.fcgi", params=params, timeout=60)
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(f"Error fetching NCBI batch at offset {start}: {e}")
                continue

            current_header, current_seq = None, []
            for line in response.text.split('\n'):
                if line.startswith('>'):
                    if current_header:
                        yield (current_header, ''.join(current_seq))
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_header:
                yield (current_header, ''.join(current_seq))

            # NCBI: max 3 req/s without an API key, 10 req/s with one.
            time.sleep(0.4 if self.api_key else 1.0)


class ViPRDownloader:
    """
    Download DENV sequences from BV-BRC (bv-brc.org).

    ViPR (viprbrc.org) was decommissioned; its data and API were folded into
    BV-BRC. The class name is kept as ViPRDownloader for backward compatibility
    with existing callers/imports.
    """

    BASE_URL = "https://www.bv-brc.org/api"
    DENV_TAXON_LINEAGE_ID = "12637"  # Dengue virus (covers all 4 serotypes)

    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({"accept": "application/json"})

    def download_dengue_genomes(self, max_genomes: int = 500) -> List[Dict]:
        """
        Fetch DENV genome metadata from BV-BRC.

        Args:
            max_genomes: Maximum number of genome records to retrieve.

        Returns:
            List of normalized genome dictionaries.
        """
        endpoint = f"{self.BASE_URL}/genome/"
        query = (
            f"eq(taxon_lineage_ids,{self.DENV_TAXON_LINEAGE_ID})"
            f"&select(genome_id,genome_name,strain,collection_date,"
            f"isolation_country,geographic_location,genome_length)"
            f"&limit({max_genomes})"
        )

        logger.info("Fetching DENV genomes from BV-BRC...")
        try:
            response = self.session.get(f"{endpoint}?{query}", timeout=60)
            response.raise_for_status()
            genomes = response.json()
            logger.info(f"Retrieved {len(genomes)} genome records from BV-BRC")

            normalized = []
            for g in genomes:
                normalized.append({
                    'genome_id': g.get('genome_id'),
                    'strain': g.get('strain') or g.get('genome_name'),
                    'collection_date': g.get('collection_date'),
                    'country': g.get('isolation_country'),
                    'location': g.get('geographic_location'),
                    'length': g.get('genome_length'),
                    'host': 'Homo sapiens',
                    'source_database': 'BV-BRC',
                })
            return normalized
        except (requests.RequestException, ValueError) as e:
            logger.error(f"Error fetching BV-BRC data: {e}")
            return []

    def download_sequence(self, genome_id: str) -> Optional[str]:
        """Download a single genome's FASTA sequence from BV-BRC by genome_id."""
        endpoint = f"{self.BASE_URL}/genome_sequence/"
        query = f"eq(genome_id,{genome_id})&http_accept=application/dna+fasta"

        try:
            response = self.session.get(f"{endpoint}?{query}", timeout=30,
                                         headers={"accept": "text/plain"})
            response.raise_for_status()
            return response.text if response.text.startswith(">") else None
        except requests.RequestException as e:
            logger.error(f"Error downloading BV-BRC sequence {genome_id}: {e}")
            return None

    def bulk_download_fasta(self, output_dir: str = "bvbrc_dengue",
                             max_genomes: int = 500) -> Generator[Tuple[str, str], None, None]:
        """
        Bulk download DENV genome sequences (as FASTA) from BV-BRC.

        Args:
            output_dir: Present for API-compatibility; not used to write files here.
            max_genomes: Cap on how many genomes to fetch.

        Yields:
            (header, sequence) tuples
        """
        genomes = self.download_dengue_genomes(max_genomes=max_genomes)
        genome_ids = [g['genome_id'] for g in genomes if g.get('genome_id')]

        if not genome_ids:
            return

        endpoint = f"{self.BASE_URL}/genome_sequence/"
        batch_size = 50

        for i in range(0, len(genome_ids), batch_size):
            batch = genome_ids[i:i + batch_size]
            id_list = ",".join(batch)
            query = f"in(genome_id,({id_list}))&http_accept=application/dna+fasta&limit({batch_size})"

            try:
                response = self.session.get(f"{endpoint}?{query}", timeout=60,
                                             headers={"accept": "text/plain"})
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(f"Error fetching BV-BRC sequence batch: {e}")
                continue

            current_header, current_seq = None, []
            for line in response.text.split('\n'):
                if line.startswith('>'):
                    if current_header:
                        yield (current_header, ''.join(current_seq))
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_header:
                yield (current_header, ''.join(current_seq))

            time.sleep(0.5)


class NextstrainDownloader:
    """
    Download curated DENV sequences and metadata from Nextstrain's public
    data bucket (data.nextstrain.org), the same files nextstrain/dengue's
    own build consumes. Files are zstd-compressed.
    """

    DATA_BASE_URL = "https://data.nextstrain.org/files/workflows/dengue"
    CHARON_URL = "https://nextstrain.org/charon/getDataset"

    def __init__(self):
        self.session = requests.Session()

    @staticmethod
    def _decompress_zst(content: bytes) -> str:
        """Decompress zstd-compressed bytes to text."""
        try:
            import zstandard
        except ImportError:
            raise RuntimeError(
                "The 'zstandard' package is required to read Nextstrain data "
                "(pip install zstandard)"
            )
        dctx = zstandard.ZstdDecompressor()
        return dctx.decompress(content, max_output_size=500_000_000).decode('utf-8', errors='replace')

    def download_serotype_data(self, serotype: int = 1) -> Dict:
        """
        Download sequences, metadata, and phylogenetic tree for a serotype.

        Args:
            serotype: DENV serotype (1-4)

        Returns:
            Dictionary with any of 'sequences' (FASTA text), 'metadata' (TSV text),
            'tree' (auspice dataset JSON) that downloaded successfully.
        """
        result = {}
        label = f"denv{serotype}"

        # Sequences (zstd-compressed FASTA)
        seq_url = f"{self.DATA_BASE_URL}/sequences_{label}.fasta.zst"
        logger.info(f"Downloading DENV-{serotype} sequences from Nextstrain...")
        try:
            response = self.session.get(seq_url, timeout=120)
            response.raise_for_status()
            result['sequences'] = self._decompress_zst(response.content)
        except (requests.RequestException, RuntimeError) as e:
            logger.error(f"Error downloading Nextstrain sequences: {e}")

        # Metadata (zstd-compressed TSV)
        meta_url = f"{self.DATA_BASE_URL}/metadata_{label}.tsv.zst"
        try:
            response = self.session.get(meta_url, timeout=120)
            response.raise_for_status()
            result['metadata'] = self._decompress_zst(response.content)
        except (requests.RequestException, RuntimeError) as e:
            logger.error(f"Error downloading Nextstrain metadata: {e}")

        # Phylogenetic tree (Auspice v2 dataset JSON, served live via charon)
        try:
            response = self.session.get(self.CHARON_URL, params={"prefix": f"/dengue/{label}"}, timeout=30)
            response.raise_for_status()
            result['tree'] = response.json()
        except requests.RequestException as e:
            logger.error(f"Error downloading Nextstrain tree: {e}")

        logger.info(f"Downloaded DENV-{serotype} data from Nextstrain")
        return result

    def parse_metadata(self, metadata_tsv: str) -> List[Dict]:
        """
        Parse Nextstrain metadata TSV into normalized genome dicts, keyed
        for later joining against sequence records by 'accession'.
        """
        genomes = []
        lines = metadata_tsv.strip().split('\n')

        if not lines:
            return genomes

        headers = lines[0].split('\t')

        for line in lines[1:]:
            fields = line.split('\t')
            if len(fields) != len(headers):
                continue

            record = dict(zip(headers, fields))
            genome = {
                'accession': record.get('accession'),
                'strain': record.get('strain'),
                'collection_date': record.get('date'),
                'country': record.get('country'),
                'region': record.get('region'),
                'genotype': record.get('genotype') or record.get('genotype_nextclade'),
                'clade': record.get('major_lineage') or record.get('minor_lineage'),
                'source_database': 'Nextstrain',
            }
            genomes.append(genome)

        return genomes


class GenBankFTPDownloader:
    """Download DENV sequences from NCBI GenBank FTP (for legacy/complete archive)."""

    BASE_URL = "ftp://ftp.ncbi.nlm.nih.gov/genbank/"

    def __init__(self):
        self.ftp = None

    def connect(self):
        """Connect to NCBI FTP server."""
        import ftplib
        try:
            self.ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
            self.ftp.login()
            logger.info("Connected to NCBI FTP server")
        except Exception as e:
            logger.error(f"Error connecting to NCBI FTP: {e}")

    def download_dengue_genbank_files(self, output_dir: str = "genbank_dengue") -> List[str]:
        """
        Download GenBank flat files containing DENV sequences.
        This is slower but provides complete historical archive.

        Args:
            output_dir: Directory to save files

        Returns:
            List of downloaded file paths
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        if not self.ftp:
            self.connect()

        logger.info("Downloading GenBank DENV files (this may take a while)...")
        downloaded_files = []

        # Note: GenBank organizes files by release (gb*.seq.gz)
        # This is a placeholder; actual implementation would iterate genbank files
        # and filter for Dengue entries

        logger.info(f"Would download to {output_path}")
        return downloaded_files


def download_all_sources(output_dir: str = "dengue_data", max_genomes: int = 500) -> Dict[str, int]:
    """
    Orchestrate downloads from all sources and write results to FASTA
    (and, for Nextstrain, TSV metadata) files under output_dir.

    Args:
        output_dir: Base directory for all downloads
        max_genomes: Cap on genomes fetched per bulk-downloadable source
                     (NCBI, BV-BRC). Nextstrain's curated files are fetched
                     in full since they're pre-filtered/curated already.

    Returns:
        Dictionary with per-source download counts.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    stats = {
        'ncbi': 0,
        'bvbrc': 0,
        'nextstrain': 0
    }

    # NCBI (E-utilities)
    ncbi_downloader = NCBIDatasetsDownloader()
    ncbi_file = output_path / "ncbi_dengue.fasta"
    count = 0
    with open(ncbi_file, 'w', encoding='utf-8') as f:
        for header, seq in ncbi_downloader.bulk_download_fasta(max_genomes=max_genomes):
            f.write(f">{header}\n{seq}\n")
            count += 1
    stats['ncbi'] = count
    logger.info(f"NCBI: {count} genomes saved to {ncbi_file}")

    # BV-BRC (successor to ViPR)
    bvbrc_downloader = ViPRDownloader()
    bvbrc_file = output_path / "bvbrc_dengue.fasta"
    count = 0
    with open(bvbrc_file, 'w', encoding='utf-8') as f:
        for header, seq in bvbrc_downloader.bulk_download_fasta(max_genomes=max_genomes):
            f.write(f">{header}\n{seq}\n")
            count += 1
    stats['bvbrc'] = count
    logger.info(f"BV-BRC: {count} genomes saved to {bvbrc_file}")

    # Nextstrain (curated, all 4 serotypes)
    nextstrain_downloader = NextstrainDownloader()
    total_nextstrain = 0
    for serotype in range(1, 5):
        ns_data = nextstrain_downloader.download_serotype_data(serotype)

        if 'sequences' in ns_data:
            seq_file = output_path / f"nextstrain_denv{serotype}.fasta"
            with open(seq_file, 'w', encoding='utf-8') as f:
                f.write(ns_data['sequences'])
            logger.info(f"Saved Nextstrain DENV-{serotype} sequences to {seq_file}")
            total_nextstrain += 1

        if 'metadata' in ns_data:
            meta_file = output_path / f"nextstrain_denv{serotype}_metadata.tsv"
            with open(meta_file, 'w', encoding='utf-8') as f:
                f.write(ns_data['metadata'])
            logger.info(f"Saved Nextstrain DENV-{serotype} metadata to {meta_file}")

    stats['nextstrain'] = total_nextstrain

    logger.info(f"Download complete: {stats}")
    return stats


if __name__ == "__main__":
    # Test downloaders
    stats = download_all_sources(output_dir="test_dengue_data", max_genomes=20)
    print(f"Download Statistics: {stats}")
