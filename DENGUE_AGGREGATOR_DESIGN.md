# Comprehensive Dengue Virus Genome Aggregator
## Architecture & Implementation Guide

**Purpose**: Unified database of all historical DENV sequences with variant annotation and outbreak tracking for Sri Lankan epidemiological context.

**Target Output**: SQLite/PostgreSQL database enabling real-time outbreak investigation by variant, lineage, geography, and temporal trends.

---

## 1. System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                      DATA SOURCE LAYER                          │
├──────────────┬──────────────┬──────────────┬──────────────┬─────┤
│ NCBI Datasets│ NCBI GenBank │   ViPR API   │  Nextstrain  │ DEN │
│   (JSON)     │   (Flat)     │   (JSON)     │   (auspice)  │ Vir │
└──────┬───────┴──────┬───────┴──────┬───────┴──────┬───────┴─────┘
       │              │              │              │
       └──────────────┴──────────────┴──────────────┘
                      │
            ┌─────────▼────────────┐
            │  PARSER/NORMALIZER   │
            │  (FASTA/GenBank/JSON)│
            └─────────┬────────────┘
                      │
       ┌──────────────┴──────────────┐
       │                             │
  ┌────▼──────────────┐   ┌──────────▼────────┐
  │  VARIANT ENGINE   │   │ PROTEIN ANNOTATOR │
  │ (SNP/Indel/Impact)│   │ (ORF translation) │
  └────┬──────────────┘   └──────────┬────────┘
       │                             │
       └──────────────┬──────────────┘
                      │
            ┌─────────▼────────────┐
            │   DATABASE LAYER     │
            │  (SQLite/PostgreSQL) │
            └─────────┬────────────┘
                      │
       ┌──────────────┴──────────────┐
       │                             │
  ┌────▼─────────────┐    ┌──────────▼────────┐
  │  QUERY INTERFACE │    │  EXPORT ENGINES   │
  │  (SQL/CLI)       │    │ (FASTA/VCF/JSON)  │
  └──────────────────┘    └───────────────────┘
```

---

## 2. Database Schema

### Core Tables

#### `genomes`
```sql
CREATE TABLE genomes (
  genome_id TEXT PRIMARY KEY,
  accession VARCHAR(20) UNIQUE NOT NULL,
  serotype INT CHECK(serotype IN (1,2,3,4)),
  lineage VARCHAR(50),
  length INT,
  sequence TEXT,
  collection_date DATE,
  collection_year INT,
  country VARCHAR(100),
  region VARCHAR(100),
  host TEXT DEFAULT 'Homo sapiens',
  source_database VARCHAR(50),  -- NCBI, ViPR, Nextstrain, etc.
  pubmed_id INT,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_serotype ON genomes(serotype);
CREATE INDEX idx_country ON genomes(country);
CREATE INDEX idx_collection_date ON genomes(collection_date);
CREATE INDEX idx_lineage ON genomes(lineage);
```

#### `variants`
```sql
CREATE TABLE variants (
  variant_id TEXT PRIMARY KEY,
  genome_id TEXT NOT NULL REFERENCES genomes(genome_id),
  position INT NOT NULL,
  ref_allele VARCHAR(100),
  alt_allele VARCHAR(100),
  variant_type VARCHAR(20),  -- SNP, INDEL, MNP
  annotation TEXT,  -- e.g., "NS1:p.A101V"
  impact VARCHAR(20),  -- HIGH, MODERATE, LOW, MODIFIER
  protein VARCHAR(50),  -- C, prM, E, NS1, NS2A, etc.
  codon_change VARCHAR(20),  -- e.g., "GCG>GTG"
  aa_change VARCHAR(20),  -- e.g., "A101V"
  is_synonymous BOOLEAN,
  is_frameshift BOOLEAN,
  is_stop_gained BOOLEAN,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_genome_variants ON variants(genome_id);
CREATE INDEX idx_position ON variants(position);
CREATE INDEX idx_protein ON variants(protein);
CREATE INDEX idx_impact ON variants(impact);
```

#### `proteins`
```sql
CREATE TABLE proteins (
  protein_id TEXT PRIMARY KEY,
  genome_id TEXT NOT NULL REFERENCES genomes(genome_id),
  protein_name VARCHAR(50),  -- C, prM, E, NS1, NS2A, NS2B, NS3, NS4A, NS4B, NS5
  start_pos INT,
  end_pos INT,
  sequence TEXT,
  function TEXT,
  pdb_id VARCHAR(10),  -- Link to 3D structure if available
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_protein_genome ON proteins(genome_id);
CREATE INDEX idx_protein_name ON proteins(protein_name);
```

#### `epitopes`
```sql
CREATE TABLE epitopes (
  epitope_id TEXT PRIMARY KEY,
  protein VARCHAR(50),
  start_pos INT,
  end_pos INT,
  epitope_sequence VARCHAR(100),
  epitope_type VARCHAR(50),  -- e.g., B-cell, T-cell CD8, T-cell CD4
  immune_response VARCHAR(100),
  source_database VARCHAR(50),  -- DENV-Ab DB, IEDB, etc.
  variants_in_epitope_count INT DEFAULT 0,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_epitope_protein ON epitopes(protein);
```

#### `lineage_assignments`
```sql
CREATE TABLE lineage_assignments (
  assignment_id TEXT PRIMARY KEY,
  genome_id TEXT NOT NULL REFERENCES genomes(genome_id),
  serotype INT,
  genotype VARCHAR(50),  -- e.g., DENV-1 Cosmopolitan, DENV-2 Asian I
  clade VARCHAR(50),  -- e.g., Clade A, Clade B
  assignment_method VARCHAR(50),  -- e.g., Nextstrain, sequence-based
  confidence FLOAT,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### `metadata`
```sql
CREATE TABLE metadata (
  metadata_id TEXT PRIMARY KEY,
  genome_id TEXT NOT NULL REFERENCES genomes(genome_id),
  key VARCHAR(100),
  value TEXT,
  source_database VARCHAR(50),
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

## 3. Data Sources & Parsing Strategy

### 3.1 NCBI Datasets API
**Endpoint**: `https://api.ncbi.nlm.nih.gov/datasets/v1/virus/taxon_tree?taxons=12637`

**Format**: JSON with hierarchical structure

**Key Fields**: 
- `accession`, `organism_name`, `virus_name`
- `collection_date`, `geo_location`
- `sequence` (or download link)
- Protein annotations (CDS features)

**Download Strategy**:
```python
GET /datasets/v1/virus/taxon_tree?taxons=12637  # All DENV
GET /datasets/v1/virus/by_accession?accesssions=NC_001474  # Specific
GET /datasets/download?taxons=12637&filename=dengue_genomes.zip  # Bulk
```

---

### 3.2 NCBI GenBank FTP
**Source**: `ftp://ftp.ncbi.nlm.nih.gov/genbank/`

**Format**: Flat-file GenBank records

**Parsing Challenge**: 
- Extract LOCUS, DEFINITION, ACCESSION, VERSION
- Parse FEATURES (CDS, misc_feature, repeat_region)
- Extract ORIGIN (sequence)

---

### 3.3 ViPR (Virus Pathogen Resource)
**Endpoint**: `http://www.viprbrc.org/brc/api/`

**Format**: JSON REST API

**Key Endpoint**:
```
GET /genome?virus_type=Flavivirus&host=Homo%20sapiens&format=json
```

**Fields**:
- `sequence_id`, `strain_name`, `collection_date`
- `country`, `segment` (for segmented viruses, N/A for DENV)
- Full sequence or download link
- Protein sequences + annotations

---

### 3.4 Nextstrain
**Endpoint**: `https://nextstrain.org/charon/v2/dengue/` (+ serotype-specific)

**Format**: Auspice JSON (phylogenetic + metadata)

**Download Strategy**:
```
GET /charon/v2/dengue/denv1/sequences.fasta
GET /charon/v2/dengue/denv1/metadata.tsv
GET /charon/v2/dengue/denv1/tree.json  # Phylogenetic context
```

**Key Data**:
- Phylogenetic tree (topology, branch lengths)
- Clade/lineage assignments
- Geographic coordinates
- Temporal constraints

---

### 3.5 DENVirDB
**Web Portal**: `https://denvirdb.isbg.ac.in/` (if available) or direct DB access

**Parsing Strategy**: Web scraping + static downloads

**Key Annotations**:
- South Asian isolates (primary focus)
- Protein physicochemical properties
- Secondary structure predictions

---

### 3.6 DENV-Ab DB & DenvInt
**DENV-Ab DB**: `https://denvabdb.bhsai.org/`
- Monoclonal antibody epitope mappings
- Neutralization profiles

**DenvInt**: PMC5648114 database
- Protein-protein interactions (dengue proteins ↔ human host)
- Mosquito-dengue interactions

---

## 4. Variant Annotation Pipeline

### 4.1 Reference Genomes
Use consensus/representative genomes for each serotype:
- **DENV-1**: NC_001474 (Hawaiian)
- **DENV-2**: NC_001438 (New Guinea C)
- **DENV-3**: NC_001469 (Sleman 78)
- **DENV-4**: NC_002640 (814669-H)

### 4.2 Alignment & Variant Calling
```
Input: Query genome (FASTA)
↓
1. Align to serotype reference (mafft --auto or minimap2)
2. Extract variants (SAMtools mpileup or bcftools)
3. Annotate impact:
   - SNP → translate codon change
   - INDEL → frameshift? stop-gained?
   - Protein annotations (NS1, E, etc.)
4. Cross-reference epitopes & known variants
↓
Output: Variant table (impact, protein, position, aa_change)
```

### 4.3 Impact Prediction Classes
- **HIGH**: frameshift, stop gained/lost, splice site
- **MODERATE**: missense (aa change), inframe indel
- **LOW**: synonymous
- **MODIFIER**: intergenic, upstream/downstream

---

## 5. Outbreak Investigation Queries

### Example Query 1: Recent DENV-2 in Sri Lanka with Novel Variants
```sql
SELECT g.accession, g.collection_date, COUNT(v.variant_id) as variant_count,
       GROUP_CONCAT(v.aa_change, ', ') as aa_changes
FROM genomes g
LEFT JOIN variants v ON g.genome_id = v.genome_id
WHERE g.serotype = 2
  AND g.country = 'Sri Lanka'
  AND g.collection_date > DATE('now', '-6 months')
  AND v.impact IN ('HIGH', 'MODERATE')
GROUP BY g.genome_id
ORDER BY g.collection_date DESC;
```

### Example Query 2: Variants in Envelope Protein (E), Last 12 Months
```sql
SELECT g.accession, g.collection_date, g.country,
       v.aa_change, v.impact, e.epitope_type
FROM genomes g
JOIN variants v ON g.genome_id = v.genome_id
LEFT JOIN epitopes e ON v.protein = e.protein 
  AND v.position BETWEEN e.start_pos AND e.end_pos
WHERE v.protein = 'E'
  AND v.impact IN ('HIGH', 'MODERATE')
  AND g.collection_date > DATE('now', '-12 months')
ORDER BY g.collection_date DESC;
```

### Example Query 3: Geographic Spread of Specific Variant (e.g., E:p.A101V)
```sql
SELECT g.country, g.region, COUNT(*) as isolates, 
       MIN(g.collection_date) as first_detected, 
       MAX(g.collection_date) as last_detected
FROM genomes g
JOIN variants v ON g.genome_id = v.genome_id
WHERE v.aa_change = 'A101V' AND v.protein = 'E'
GROUP BY g.country, g.region
ORDER BY COUNT(*) DESC;
```

---

## 6. Implementation Roadmap

### Phase 1: Foundation (Week 1)
- [ ] Database schema setup (SQLite initially, migrate to PostgreSQL if needed)
- [ ] NCBI Datasets API downloader + parser
- [ ] Basic genome ingestion pipeline
- [ ] Serotype detection (sequence signature or metadata parsing)

### Phase 2: Multi-Source Integration (Week 2-3)
- [ ] ViPR API integration
- [ ] Nextstrain auspice JSON parser
- [ ] DENVirDB web scraping or direct import
- [ ] Deduplication logic (same genome from multiple sources)

### Phase 3: Annotation Engine (Week 3-4)
- [ ] Reference genome setup (4 DENV serotypes)
- [ ] Alignment pipeline (mafft + variant calling)
- [ ] Protein translation + impact prediction
- [ ] Epitope cross-reference (DENV-Ab DB, DenvInt)

### Phase 4: Outbreak Tools (Week 4-5)
- [ ] Query CLI interface
- [ ] Export engines (FASTA, VCF, JSON)
- [ ] Temporal clustering analysis
- [ ] Geographic heat map generation (optional)

### Phase 5: Validation & Deployment (Week 5-6)
- [ ] Sanity checks (genome counts, variant distributions)
- [ ] Comparison with published lineages (Nextstrain, ViPR)
- [ ] Documentation + example queries
- [ ] Local isolate integration pathway (for future use)

---

## 7. Software Dependencies

```
Core:
- Python 3.9+
- sqlite3 (built-in) or psycopg2 (PostgreSQL)
- requests (HTTP downloads)
- pandas (data wrangling)
- biopython (FASTA/GenBank parsing)

Alignment & Annotation:
- mafft (sequence alignment)
- samtools (BAM manipulation)
- bcftools (variant calling)
- minimap2 (fast alignment, optional)

Optional (for structure prediction):
- esm (ESMFold, protein 3D structures)
- foldseek (structure search, relevant to your polymerase manifold idea)

Visualization (post-processing):
- plotly / matplotlib (interactive plots)
- ete3 (phylogenetic trees)
```

---

## 8. Outputs

### 8.1 Exportable Formats
- **FASTA**: Filtered by serotype, lineage, country, date range
- **VCF**: Variant calls with annotations (impact, protein, epitope overlap)
- **GFF**: Genomic feature coordinates (proteins, epitopes, variants)
- **JSON**: Structured metadata + variants (web-ready)
- **TSV**: Query results for epidemiological analysis

### 8.2 Standard Query Set
- "Genomes by country & date range"
- "High-impact variants in NS1 protein (last 6 months)"
- "Geographic spread of variant X"
- "Lineage-specific variants (e.g., DENV-2 Asian I clade A)"
- "Epitope-overlapping variants (immune escape risk)"

---

## 9. Integration with Your NextSeq Workflow

Once complete, this aggregator enables:
1. **Reference genome selection**: Align your patient isolates to matching serotype/lineage reference
2. **Variant comparison**: Call SNPs/indels, compare to aggregated database
3. **Outbreak context**: Place your isolates in global phylogeny + identify locally-novel variants
4. **Structural impact**: Predict protein consequences (feed to ESMFold for 3D structure)
5. **Publication-ready figures**: Phylogenetic trees + variant timelines

---

## 10. Notes for Sri Lankan Context

- Dengue historically endemic (DENV-1, -2 most common; DENV-3, -4 occasional)
- Recent outbreaks: watch for **recombinant lineages** (e.g., DENV-2 Asian I × American genotypes)
- **Geographic clustering**: Colombo/Western Region likely dominates sequence data; ensure rural/endemic zones represented
- **Temporal signal**: Collection date metadata critical for inferring transmission chains
- **Host factors**: Immunocompromised (cancer, transplant) patients may shed variants longer—flag if available

---

## References & APIs

| Source | API/URL | Format | Rate Limit |
|--------|---------|--------|-----------|
| NCBI Datasets | api.ncbi.nlm.nih.gov | JSON | 10 req/sec |
| NCBI GenBank FTP | ftp.ncbi.nlm.nih.gov | Flat file | N/A (FTP) |
| ViPR | viprbrc.org/brc/api | JSON | 1 req/sec (check docs) |
| Nextstrain | nextstrain.org/charon | JSON/FASTA | Generous |
| DENVirDB | denvirdb.isbg.ac.in | Web (scrape) | N/A |
| DENV-Ab DB | denvabdb.bhsai.org | Web (scrape) | N/A |

