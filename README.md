# miRNA Experiment Search

A flexible Python pipeline to search public NCBI metadata sources for experiments involving selected miRNAs, cell lines, species, assay types, perturbation types, and control conditions.

The tool currently searches **GEO** and **SRA** through the NCBI Entrez E-utilities API and produces clean candidate tables plus summary statistics.

---

## Overview

This repository is designed for exploratory discovery of public experiments related to miRNAs. It supports searches such as:

- miRNA perturbation RNA-seq experiments
- AGO/CLIP/CLASH studies for miRNA-target interactions
- qPCR or reporter-assay validation studies
- experiments in specific cell lines
- studies restricted to a species
- studies with control-like metadata

The pipeline uses user-friendly canonical filter names, such as `RNA-seq`, `CLIP`, `CLASH`, `qPCR`, `overexpression`, or `knockdown`, and expands them internally into broader keyword dictionaries.

---

## Key Features

- Search by miRNA list, cell-line list, or both
- Search GEO and/or SRA
- Filter by species, for example `Homo sapiens`
- Filter by experiment type:
  - `RNA-seq`
  - `small_RNA-seq`
  - `CLIP`
  - `CLASH`
  - `qPCR`
  - `microarray`
  - `proteomics`
  - `reporter_assay`
- Filter by perturbation type:
  - `overexpression`
  - `knockdown`
  - `knockout`
  - `perturbation_any`
- Optional control keyword requirement
- Optional date filtering
- Built-in exclusion of precursor/primary miRNA-focused studies
- Generates:
  - candidate experiment table
  - summary statistics
  - summary by miRNA
  - summary by accession

---

## Repository Structure

```text
mirna-experiment-search/
├── config/
│   └── search.yaml
├── data/
│   ├── mirnas.txt
│   └── cell_lines.txt
├── results/
├── src/
│   ├── annotate.py
│   ├── entrez.py
│   ├── io_utils.py
│   ├── query_builder.py
│   ├── summarize.py
│   └── vocab.py
├── run_search.py
├── pyproject.toml
├── uv.lock
└── README.md
```

---

## Installation

This project uses [`uv`](https://docs.astral.sh/uv/) for environment and dependency management.

From inside the cloned repository:

```bash
uv sync
```

To check that the environment works:

```bash
uv run python --version
```

---

## Input Files

### miRNA list

Provide a plain text file with one miRNA per line.

Example: `data/mirnas.txt`

```text
hsa-miR-22-3p
hsa-miR-192-5p
hsa-miR-200c-3p
```

### Cell-line list

Provide a plain text file with one cell line per line.

Example: `data/cell_lines.txt`

```text
HEK293T
HEK293
HeLa
```

Both files are optional. The pipeline can search by miRNA only, cell line only, experiment type only, or any combination.

---

## Configuration

The main configuration file is:

```text
config/search.yaml
```

Example:

```yaml
mirnas_file: "data/mirnas.txt"
cell_lines_file: "data/cell_lines.txt"

sources:
  - GEO
  - SRA

species: "Homo sapiens"

experiment_types:
  - RNA-seq

perturbation_types: []

require_control: false

date_from: null
date_to: null

retmax: 20
sleep_seconds: 0.34

exclude_precursor_terms: true

outdir: "results"
```

The config file contains only high-level user choices. Synonyms and related keywords are defined internally in `src/vocab.py`.

For example:

```yaml
experiment_types:
  - RNA-seq
```

is expanded internally to terms such as:

```text
RNA-seq
RNA seq
RNA sequencing
mRNA-seq
transcriptome
transcriptomic
```

---

## Running the Pipeline

Run with the default config:

```bash
uv run python run_search.py --config config/search.yaml
```

Run with command-line overrides:

```bash
uv run python run_search.py \
  --config config/search.yaml \
  --sources GEO \
  --species "Homo sapiens" \
  --experiment-types RNA-seq CLIP \
  --perturbation-types overexpression knockdown \
  --retmax 10 \
  --outdir results/test_run
```

On Windows PowerShell:

```powershell
uv run python run_search.py `
  --config config/search.yaml `
  --sources GEO `
  --species "Homo sapiens" `
  --experiment-types RNA-seq CLIP `
  --perturbation-types overexpression knockdown `
  --retmax 10 `
  --outdir results/test_run
```

---

## Example Searches

### Search human RNA-seq studies for selected miRNAs

```bash
uv run python run_search.py \
  --mirnas-file data/mirnas.txt \
  --species "Homo sapiens" \
  --experiment-types RNA-seq
```

### Search HEK-related CLIP and CLASH experiments

```bash
uv run python run_search.py \
  --cell-lines-file data/cell_lines.txt \
  --species "Homo sapiens" \
  --experiment-types CLIP CLASH
```

### Search miRNA overexpression experiments

```bash
uv run python run_search.py \
  --mirnas-file data/mirnas.txt \
  --species "Homo sapiens" \
  --perturbation-types overexpression
```

### Require control-like metadata

```bash
uv run python run_search.py \
  --mirnas-file data/mirnas.txt \
  --species "Homo sapiens" \
  --experiment-types RNA-seq \
  --require-control
```

### Search with date bounds

```bash
uv run python run_search.py \
  --mirnas-file data/mirnas.txt \
  --species "Homo sapiens" \
  --experiment-types RNA-seq \
  --date-from 2018 \
  --date-to 2026
```

---

## Outputs

The pipeline writes output files to the configured output directory.

Default:

```text
results/
```

### `candidate_experiments.tsv`

Main output table with one row per candidate record.

Typical columns include:

| Column | Description |
|---|---|
| `query_mirna` | miRNA used in the query |
| `query_cell_line` | cell line used in the query |
| `source` | GEO or SRA |
| `uid` | Entrez UID |
| `accession` | GEO/SRA accession when available |
| `summary` | metadata summary, if available |
| `species_or_organism` | organism metadata |
| `matched_experiment_types` | inferred experiment categories |
| `matched_perturbation_types` | inferred perturbation categories |
| `has_perturbation_keyword` | whether perturbation terms were found |
| `has_overexpression_keyword` | whether overexpression terms were found |
| `has_knockdown_keyword` | whether knockdown terms were found |
| `has_knockout_keyword` | whether knockout terms were found |
| `has_control_keyword` | whether control-like terms were found |
| `inferred_cell_lines` | cell lines inferred from metadata |
| `url` | NCBI URL |
| `title` | record title |

The internal Entrez query is not included in the final candidate table.

### `summary_statistics.tsv`

Compact count table with overall statistics, such as:

- total candidate records
- unique accessions
- records by source
- records with controls
- records by experiment type
- records by perturbation type
- records by inferred cell line

### `summary_by_mirna.tsv`

One row per queried miRNA with:

- number of records
- unique accessions
- GEO/SRA counts
- control counts
- inferred cell lines
- matched experiment types
- matched perturbation types

### `summary_by_accession.tsv`

One row per accession with:

- source
- accession
- number of matched rows
- matched miRNAs
- matched cell lines
- inferred experiment and perturbation types
- title
- URL

---

## Entrez Search Implementation

Searches are performed using the NCBI Entrez Programming Utilities, specifically:

- `esearch` to retrieve matching record UIDs
- `esummary` to retrieve metadata for those UIDs

For each combination of miRNA, cell line, source, and filters, the pipeline builds an Entrez query using:

- miRNA aliases
- optional cell-line terms
- species terms
- experiment-type synonyms
- perturbation-type synonyms
- optional control terms
- optional date range
- optional precursor/primary miRNA exclusion terms

The query is submitted to either:

- GEO via Entrez database `gds`
- SRA via Entrez database `sra`

A short delay is applied between requests to avoid excessive request rates.

---

## Built-in Vocabulary

The controlled vocabulary is defined in:

```text
src/vocab.py
```

This includes dictionaries for:

- experiment types
- perturbation types
- control terms
- species terms
- precursor-exclusion terms
- common cell-line aliases

This design keeps the user config simple while allowing the code to expand high-level categories into comprehensive search terms.

---

## Limitations

This tool searches metadata, not raw experimental data. Therefore:

- GEO/SRA metadata can be incomplete or noisy
- cell-line inference is keyword-based and approximate
- control detection is keyword-based and approximate
- a candidate hit does not guarantee the experiment is directly usable
- manual curation is recommended before downstream biological analysis

---

## Recommended Workflow

1. Prepare `mirnas.txt` and/or `cell_lines.txt`
2. Select filters in `config/search.yaml`
3. Run the pipeline
4. Inspect `candidate_experiments.tsv`
5. Use `summary_statistics.tsv` and `summary_by_mirna.tsv` for triage
6. Manually verify top candidate studies in GEO/SRA

---

## Development

Run the pipeline:

```bash
uv run python run_search.py --config config/search.yaml
```

Check Git status:

```bash
git status
```

Commit changes:

```bash
git add .
git commit -m "Describe your change"
```

---

## Citation and Usage

This repository is intended as a metadata search and triage tool for miRNA-related public experiments. It should be used as a first-pass discovery workflow before manual curation and downstream analysis.
