# GEO miRNA OE/KO RNA-seq Search

Search GEOmetadb for human RNA-seq GEO series where mature miRNAs are
overexpressed or knocked out.

## Install

```bash
cd /home/dtzim01/geo_rna_perturbation_search
conda env create -f environment.yml
conda activate geo-rna-perturbation
```

If the environment already exists, only run:

```bash
cd /home/dtzim01/geo_rna_perturbation_search
conda activate geo-rna-perturbation
```

The script uses `GEOmetadb.sqlite` in this directory. If the file is missing,
it downloads it automatically. If it exists, it reuses it.

## miRNA List

Input is a plain text file with one mature miRNA per line:

```text
hsa-miR-21-5p
hsa-let-7a-5p
```

Blank lines are ignored. Text after `#` is treated as a comment.

## Run All Dates

Search overexpression and knockout:

```bash
Rscript search_geo_rna_perturbations.R \
  --targets-file mature_mirnas_test.txt \
  --mode both \
  --cpus 8 \
  --output mature_mirna_oe_ko_geo.tsv \
  --ambiguous-output mature_mirna_oe_ko_ambiguous.tsv
```

## Output

Main output: clear human GEO series-level hits, one row per mature miRNA + GSE.

Ambiguous output: mixed/ambiguous rows excluded from the main output.

## Flags

- `--targets-file`: required text file with one mature miRNA per line.
- `--mode`: `both`, `oe`, or `ko`. Default: `both`.
- `--date`: optional comma-separated GEO submission years, e.g. `2023,2024`.
- `--cpus`: number of parallel workers. For ~600 miRNAs, start with `4` or `8`.
- `--output`: main clean TSV.
- `--ambiguous-output`: optional TSV for ambiguous rows.
