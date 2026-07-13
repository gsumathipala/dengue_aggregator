# Dengue Virus Genome Aggregator

**A comprehensive bioinformatics pipeline to download, parse, and annotate all historical DENV sequences from multiple authoritative sources into a unified, queryable database for outbreak investigation.**

---

## Overview

This aggregator consolidates dengue virus genome sequences from:
- **NCBI Datasets API** (official NCBI GenBank records)
- **Virus Pathogen Resource (ViPR)** (curated viral genomes)
- **Nextstrain** (real-time phylogenetic tracking)
- **DENVirDB** (South Asian isolate focus)
- **DENV-Ab DB** (epitope mapping)
- **DenvInt** (host-virus interactions)

It produces a unified SQLite database enabling:
1. **Rapid queries** by serotype, country, lineage, collection date
2. **Variant tracking** (SNPs, indels, functional impact)
3. **Epidemic context** (geographic spread, temporal clustering)
4. **Integration with local isolates** (from your NextSeq runs)

---

## Installation

### Prerequisites
- Python 3.9+
- SQLite3 (built-in)
- Optional: mafft, minimap2, samtools (for variant annotation)

### Setup

```bash
# Clone or download the aggregator
cd dengue_aggregator

# Install Python dependencies
pip install requests pandas biopython

# Optional: For advanced sequence alignment
pip install mafft minimap2

# Initialize database
python main.py --db dengue_genomes.db
```

---

## Quick Start

### 1. Full Pipeline (Download + Parse + Annotate)

```bash
python main.py \
  --db dengue_genomes.db \
  --data-dir dengue_data \
  -v  # verbose logging
```

**What happens:**
- Downloads from NCBI Datasets, Nextstrain, ViPR
- Parses FASTA/GenBank/JSON formats
- Quality control checks (sequence length, GC content, etc.)
- Inserts into SQLite database
- (Optional) Calls variants and predicts impact

**Expected output:**
```
[PHASE 1] Downloading from all data sources...
  NCBI Datasets: 5,234 genomes
  Nextstrain DENV-1: 1,892 sequences
  Nextstrain DENV-2: 2,104 sequences
  ...

[PHASE 2] Parsing downloaded files...
  Parsing nextstrain_denv1.fasta...
    ✓ Inserted 1,234 genomes from nextstrain_denv1.fasta
  ...

[PHASE 3] Annotating variants...
  (Requires reference genomes - can be added later)

========================================
EXECUTION SUMMARY
========================================
Downloaded:        10,000 sequences
Parsed:             9,876 sequences
Passed QC:          9,654 sequences
Failed QC:            222 sequences
Inserted to DB:     9,654 sequences
Errors:               0 issues
========================================
```

### 2. Skip Download (Use Existing Files)

```bash
python main.py \
  --skip-download \
  --db dengue_genomes.db \
  --data-dir dengue_data
```

### 3. Query by Serotype

```bash
python main.py \
  --db dengue_genomes.db \
  --query-serotype 2 \
  --query-country "Sri Lanka"
```

---

## Usage Examples

### Command Line

#### Example 1: Update DENV-2 Reference Data
```bash
python main.py \
  --db dengue_v2.db \
  --data-dir nextstrain_v2_only \
  --skip-annotate \
  -v
```

#### Example 2: Query Recent Isolates
```bash
python main.py \
  --db dengue_genomes.db \
  --query-country "Sri Lanka"
```

---

## Python API Usage

### Direct Database Access

```python
from db_schema import DengueDB
from parsers import FASTAParser, QualityControl

# Initialize database
db = DengueDB("dengue_genomes.db")

# Query by serotype
genomes_denv2 = db.query_genomes({'serotype': 2})

# Query by location and date
sri_lanka_recent = db.query_genomes({
    'country': 'Sri Lanka',
    'after_date': '2023-01-01',
    'before_date': '2024-12-31'
})

for genome in sri_lanka_recent:
    print(f"{genome['accession']}: {genome['collection_date']} from {genome['country']}")

db.close()
```

### Parsing Local FASTA Files

```python
from parsers import FASTAParser, QualityControl, GenomeRecord

# Parse your NextSeq/local sequencing results
for genome in FASTAParser.parse_file("my_samples.fasta"):
    # Check quality
    passed, issues = QualityControl.validate_genome(genome)
    
    if passed:
        print(f"✓ {genome.accession}: Valid DENV sequence")
    else:
        print(f"✗ {genome.accession}: {issues}")
```

### Variant Annotation

```python
from variant_annotator import VariantAnnotator, DENV_REFERENCES

# Get reference sequence
ref_seq = # ... download from NCBI NC_001474 (DENV-1)

# Initialize annotator
annotator = VariantAnnotator(ref_seq, serotype=1)

# Call variants
query_seq = "AGTTGTTAGGA..."  # Your query sequence
variants = annotator.call_variants(query_seq)

for variant in variants:
    print(f"{variant.protein}: {variant.aa_change} ({variant.impact})")
```

---

## Database Schema

### Core Tables

#### `genomes`
Primary table with all genome metadata:
```sql
SELECT accession, serotype, lineage, country, collection_date, source_database
FROM genomes
WHERE serotype = 2 AND country = 'Sri Lanka'
ORDER BY collection_date DESC;
```

#### `variants`
Variant calls with functional annotations:
```sql
SELECT g.accession, v.protein, v.aa_change, v.impact
FROM genomes g
JOIN variants v ON g.genome_id = v.genome_id
WHERE v.impact = 'HIGH' AND g.serotype = 2
ORDER BY v.position;
```

#### `proteins`
Annotated protein sequences (C, prM, E, NS1-NS5):
```sql
SELECT protein_name, COUNT(*) as genes
FROM proteins
WHERE genome_id = 'DENV_xxx'
GROUP BY protein_name;
```

#### `epitopes`
B-cell and T-cell epitope mappings (from DENV-Ab DB):
```sql
SELECT * FROM epitopes
WHERE protein = 'E' AND epitope_type = 'B-cell neutralizing';
```

---

## Example Queries for Outbreak Investigation

### Query 1: Recent High-Impact Variants (Last 6 Months, DENV-2)
```sql
SELECT 
    g.accession, 
    g.collection_date, 
    g.country, 
    COUNT(v.variant_id) as impact_variants,
    GROUP_CONCAT(v.aa_change, '; ') as changes
FROM genomes g
LEFT JOIN variants v ON g.genome_id = v.genome_id
WHERE g.serotype = 2 
    AND g.collection_date > date('now', '-6 months')
    AND v.impact IN ('HIGH', 'MODERATE')
GROUP BY g.genome_id
ORDER BY g.collection_date DESC;
```

### Query 2: Geographic Spread of Lineage-Specific Variants
```sql
SELECT 
    g.country, 
    COUNT(DISTINCT g.genome_id) as isolates,
    MIN(g.collection_date) as first_detected,
    MAX(g.collection_date) as last_detected
FROM genomes g
JOIN variants v ON g.genome_id = v.genome_id
WHERE g.lineage = 'Asian I' 
    AND v.aa_change = 'E:A101V'
GROUP BY g.country
ORDER BY COUNT(*) DESC;
```

### Query 3: Variants in Envelope Protein (Immune Escape Risk)
```sql
SELECT 
    g.accession,
    g.collection_date,
    g.country,
    v.aa_change,
    e.epitope_type as epitope_overlap
FROM genomes g
JOIN variants v ON g.genome_id = v.genome_id
LEFT JOIN epitopes e ON v.protein = e.protein 
    AND v.protein_position BETWEEN e.start_pos AND e.end_pos
WHERE v.protein = 'E' 
    AND g.collection_date > date('now', '-12 months')
ORDER BY g.collection_date DESC;
```

---

## Data Sources & Attribution

Sources actually integrated in this pipeline:

| Source | Coverage | Strength | Contact |
|--------|----------|----------|---------|
| NCBI Datasets/GenBank | Global, curated | Official GenBank records | ncbi.nlm.nih.gov |
| Nextstrain | Real-time, phylo | Phylogenetic context, per-serotype metadata | nextstrain.org |
| BV-BRC (Bacterial and Viral Bioinformatics Resource Center) | Curated annotations | NIAID-funded resource, successor to ViPR/PATRIC | bv-brc.org |
| IEDB (Immune Epitope Database) | B-/T-cell epitope mappings | NIAID-funded resource | iedb.org |

All four are NIH/NIAID-funded or INSDC public sequence archives intended for
open reuse; sequence and epitope data pulled from them is not subject to
additional restriction beyond attribution. **GISAID is not used anywhere in
this pipeline** and should not be added to a repository that redistributes
data — GISAID's terms of use prohibit redistributing its raw records outside
GISAID itself.

Note: this repository's own database (`dengue_genomes.db`) and downloaded raw
files (`dengue_data/`) are intentionally **not** committed to version control
(see `.gitignore`) — they're large, fully regenerable derived data. Rebuild
them locally with `main.py` / `annotate_variants.py` / `build_variant_catalog.py`
/ `load_epitopes.py` as described above.

---

## Integration with Your Workflows

### NextSeq 1000 Integration (Future Phase)
```python
# After your NGS run, align to aggregator references
from variant_annotator import VariantAnnotator

# Load your patient isolate
patient_fasta = "Sample_R1_001_trimmed.fasta"

# Get reference from aggregator
ref_genome = db.query_genomes({'serotype': 2, 'lineage': 'Asian I'})[0]

# Annotate
annotator = VariantAnnotator(ref_genome['sequence'], serotype=2)
patient_variants = annotator.call_variants(patient_fasta)

# Compare to epidemiological context
print(f"Your isolate carries {len(patient_variants)} variants")
print(f"Of these, {sum(1 for v in patient_variants if v.impact == 'HIGH')} are HIGH impact")
```

### Phylogenetic Analysis (Downstream)
```bash
# Export FASTA for phylogenetics (e.g., IQ-TREE)
sqlite3 dengue_genomes.db "SELECT '>'||accession || E'\n' || sequence FROM genomes WHERE serotype=2 AND collection_date > '2023-01-01';" > denv2_recent.fasta

# Build tree
iqtree -s denv2_recent.fasta -m GTR+G -bb 1000

# Visualize in iTOL or Figtree
```

---

## Troubleshooting

### Issue: "No FASTA files found"
**Cause**: Download phase didn't complete or files saved to wrong directory.
**Solution**: Check `--data-dir`, verify internet connection, re-run with `--skip-parse`.

### Issue: "QC failed for many genomes"
**Cause**: Downloaded sequences may include incomplete/partial genomes.
**Solution**: This is normal; QC filters out short/contaminated sequences. Check `DENV_EXPECTED_LENGTH_RANGE` in `quality_control.py`.

### Issue: "TypeError: 'NoneType' object is not iterable"
**Cause**: Parsing failed for a specific format.
**Solution**: Check `--verbose` output, inspect the problematic file manually, file an issue.

---

## Performance Notes

| Dataset Size | Download Time | Parse Time | Database Size |
|--------------|---------------|------------|---------------|
| 5,000 seqs | ~5 min | ~2 min | ~150 MB |
| 10,000 seqs | ~10 min | ~5 min | ~300 MB |
| 50,000 seqs+ | ~30-60 min | ~15-30 min | ~1.5 GB |

**Tips for large datasets:**
- Use `--skip-download` to avoid re-fetching
- Filter by date/country in download params (if supported by source)
- Run on a machine with >8 GB RAM
- Consider PostgreSQL for >100k genomes

---

## Contributing & Improvements

### Roadmap
- [ ] PostgreSQL backend for large-scale deployments
- [ ] Web interface (Flask/React) for queries
- [ ] Integration with GISAID (restricted access)
- [ ] Real-time email alerts for novel variants
- [ ] Structure prediction (ESMFold) for E protein
- [ ] Recombination detection (PANGO lineage mapping)

### Code Structure
```
dengue_aggregator/
├── main.py                  # CLI entry point
├── db_schema.py             # Database initialization
├── downloaders.py           # Source-specific downloaders
├── parsers.py               # Format parsers (FASTA, GenBank, JSON)
├── variant_annotator.py     # Variant calling + impact prediction
├── DESIGN.md                # Detailed architecture
└── README.md                # This file
```

---

## References

1. **NCBI Dengue Taxonomy**: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/12637/
2. **Nextstrain Dengue**: https://nextstrain.org/dengue
3. **ViPR Flavivirus**: http://www.viprbrc.org/brc/home.spg?decorator=vipr
4. **DENVirDB**: https://journals.lww.com/jvbd/fulltext/2014/51020/denvirdb__a_web_portal_of_dengue_virus_sequence.2.aspx
5. **DENV-Ab DB**: https://denvabdb.bhsai.org/
6. **DenvInt**: PMC5648114 (doi: 10.1186/s13015-017-0115-0)

---

## License & Citation

The code in this repository is licensed under the [MIT License](LICENSE).

That covers the pipeline code only. The underlying sequence and epitope data
is sourced from NCBI/GenBank, Nextstrain, BV-BRC, and IEDB (see **Data
Sources & Attribution** above) and remains subject to those sources' own
terms; this project does not claim any additional rights over it.

This aggregator is provided as-is for research and outbreak investigation.

**Citation:**
```bibtex
@software{dengue_aggregator_2024,
  author = {gsumathipala},
  title = {Comprehensive Dengue Virus Genome Aggregator},
  year = {2024},
  url = {https://github.com/gsumathipala/dengue_aggregator}
}
```

---

## Contact & Support

For questions, issues, or feature requests:
- File an issue on GitHub
- Contact your virology team lead
- Reach out to NCI LSUS bioinformatics core

---

**Version 1.0** | Last Updated: July 2026
