# Dengue Sequence Analyzer (Web App)

A local web interface built on top of the aggregator: paste or upload one DENV
sequence (complete genome or partial fragment) and get back serotype
identification, per-variant annotation, known-variant/epitope
cross-referencing, and nearest-neighbor placement against the ~12,000-genome
`dengue_genomes.db` database, rendered as an illustrated report.

It's the read-only, single-sequence counterpart to the batch ingestion
pipeline (`main.py` / `annotate_variants.py`) — nothing you submit here is
written back to the database.

---

## Setup

Requires the aggregator database to already exist and be annotated (see the
main `README.md`: `main.py`, then `annotate_variants.py`, then
`build_variant_catalog.py` and `load_epitopes.py`).

```bash
# Install the extra dependency (rapidfuzz/numpy are already required by the
# rest of the pipeline; see requirements-webapp.txt)
pip install flask

# One-time (and again any time the `variants` table changes): build the
# nearest-neighbor SNP fingerprint cache
python build_fingerprints.py --db dengue_genomes.db

# Run
python app.py
```

Then open **http://127.0.0.1:5000/**. Binds to localhost only, no auth — this
is a local research tool, not something to expose beyond your own machine/LAN.

---

## Using it

1. Paste a FASTA record (with or without a `>header`) or raw nucleotide text
   into the box, upload a `.fasta`/`.fa`/`.txt` file, or click **Load
   example** to pull a random real genome of a given serotype straight out of
   the database.
2. Click **Analyze**.
3. If your input contains more than one FASTA record, only the first is
   analyzed — a warning banner names it and says how many others were
   ignored.

### What you get back

- **Serotype call** — with a confidence bar chart showing k-mer containment
  against all 4 RefSeq references, so you can see the evidence, not just the
  winner.
- **Coverage** — how much of the reference genome your sequence actually
  spans (important for partial/fragment input).
- **Genome map** — an SVG track: real protein-coordinate blocks (C, prM, E,
  NS1, NS2A, NS2B, NS3, NS4A, 2K, NS4B, NS5), coverage shading, variant
  markers colored by impact, and known-epitope regions. Hover for detail.
- **Variant table** — every called variant with protein/codon/AA change,
  impact, and, for anything already seen in the database: population
  frequency, first/last observed date, countries, and lineages. Anything not
  found in `variant_catalog` is flagged **novel**.
- **Nearest neighbors** — the closest genomes in the database by shared SNP
  profile, with country/date/lineage context.
- CSV (variant table) and JSON (full report) export, generated client-side.

---

## How it works (for the curious / for maintenance)

| Step | How |
|---|---|
| Serotype + strand orientation | Exact k-mer (k=21) containment against each of the 4 serotype references, tried in both orientations. The *margin* over the runner-up serotype is the real signal — absolute containment naturally runs low (~15-25%) even for a correct call, since serotypes diverge only ~5-8% from their own RefSeq but a single mismatch kills a whole k-mer. |
| Localization | Before alignment, the query is seed-matched to estimate where in the reference it starts. This matters because naive global edit distance between a short fragment and the full-length reference degenerates into scattering the query across non-contiguous matching bases (cheap with a 4-letter alphabet) instead of aligning it to the one contiguous region it came from. See `align_query()` in `sequence_analysis.py`. |
| Variant calling | `variant_annotator.VariantAnnotator.call_variants()` — the same engine the batch pipeline uses, given the localized reference window. |
| Known-variant lookup | In-memory `variant_catalog` cache, keyed by (serotype, position, ref, alt). |
| Epitope overlap | In-memory `epitopes` cache, via the existing `annotate_with_epitopes()`. |
| Nearest neighbors | Precomputed per-serotype SNP fingerprint matrices (`fingerprints.py`) compared with one vectorized numpy op. A live query against the raw `variants` table (8M+ rows) was benchmarked at ~12s per request — too slow for interactive use — hence the precomputed cache. |

## Files

| File | Purpose |
|---|---|
| `app.py` | Flask entry point: `/`, `/api/analyze`, `/api/example` |
| `sequence_analysis.py` | The analysis pipeline (`analyze_sequence()`) |
| `fingerprints.py` | Nearest-neighbor SNP fingerprint cache: build/load/query |
| `build_fingerprints.py` | CLI to (re)build that cache |
| `templates/index.html`, `static/app.js`, `static/style.css` | The page |
| `requirements-webapp.txt` | Just `flask` |

## Caveats

- **Ephemeral by design** — no submitted sequence or result is ever written
  to disk or to `dengue_genomes.db` (every connection is opened with
  `PRAGMA query_only`).
- **Not for clinical/diagnostic use.** It's a research/exploration tool.
- Very short (<~21bp) or highly divergent fragments may fail to localize
  against the reference; the tool falls back to whole-genome alignment and
  flags this with a warning rather than failing silently.
- Rerun `build_fingerprints.py` after re-running `annotate_variants.py` or
  otherwise changing the `variants` table — the nearest-neighbor cache goes
  stale otherwise.

## Troubleshooting

- **"Fingerprint cache not found" at startup** — run `python
  build_fingerprints.py` first.
- **`ModuleNotFoundError`** for flask/rapidfuzz/numpy — make sure you're
  running `app.py` with the same Python environment the rest of the pipeline
  uses (the one with `rapidfuzz` already installed).
- **"Input is larger than X bytes"** — you likely pasted/uploaded a
  multi-genome batch FASTA rather than a single sequence; only the first
  record is ever analyzed, so trim the file down to the one you want.
