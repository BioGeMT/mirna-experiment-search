#!/usr/bin/env python3

from pathlib import Path
from itertools import product

import pandas as pd

from src.annotate import normalize_records
from src.entrez import search_entrez
from src.io_utils import load_yaml, read_lines_file
from src.query_builder import build_query
from src.summarize import (
    make_summary_by_accession,
    make_summary_by_mirna,
    make_summary_statistics,
)
from src.vocab import CONTROL_TERMS


def as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def search_combination(
    source,
    mirna,
    cell_line,
    config,
):
    query = build_query(
        mirna=mirna,
        cell_line=cell_line,
        species=config.get("species"),
        experiment_types=as_list(config.get("experiment_types")),
        perturbation_types=as_list(config.get("perturbation_types")),
        require_control=bool(config.get("require_control", False)),
        control_terms=CONTROL_TERMS,
        exclude_precursor_terms=bool(config.get("exclude_precursor_terms", True)),
        date_from=config.get("date_from"),
        date_to=config.get("date_to"),
    )

    if not query:
        raise ValueError(
            "Empty query generated. Provide at least one search condition "
            "(miRNA, cell line, species, experiment type, or perturbation type)."
        )

    print(f"\nSearching {source}")
    print(f"  miRNA: {mirna or '-'}")
    print(f"  cell line: {cell_line or '-'}")
    print(f"  query: {query}")

    records = search_entrez(
        source=source,
        query=query,
        retmax=int(config.get("retmax", 20)),
        sleep_seconds=float(config.get("sleep_seconds", 0.34)),
    )

    rows = normalize_records(
        records,
        query_mirna=mirna,
        query_cell_line=cell_line,
    )

    return rows


def apply_post_filters(df, config):
    """
    Post-filter annotated rows.

    Query-time filtering is broad. These filters apply after metadata annotation.
    """
    if df.empty:
        return df

    out = df.copy()

    if config.get("require_control", False):
        out = out[out["has_control_keyword"] == True]

    selected_experiments = set(as_list(config.get("experiment_types")))
    if selected_experiments:
        out = out[
            out["matched_experiment_types"]
            .fillna("")
            .apply(
                lambda x: bool(
                    selected_experiments.intersection(
                        set(v for v in str(x).split(";") if v)
                    )
                )
            )
        ]

    selected_perturbations = set(as_list(config.get("perturbation_types")))
    if selected_perturbations:
        out = out[
            out["matched_perturbation_types"]
            .fillna("")
            .apply(
                lambda x: bool(
                    selected_perturbations.intersection(
                        set(v for v in str(x).split(";") if v)
                    )
                )
            )
        ]

    return out


def main():
    config_path = "config/search.yaml"
    config = load_yaml(config_path)

    outdir = Path(config.get("outdir", "results"))
    outdir.mkdir(parents=True, exist_ok=True)

    mirnas = read_lines_file(config.get("mirnas_file"))
    cell_lines = read_lines_file(config.get("cell_lines_file"))

    sources = as_list(config.get("sources")) or ["GEO", "SRA"]

    # Allow searches without miRNA or without cell line.
    mirna_values = mirnas if mirnas else [None]
    cell_line_values = cell_lines if cell_lines else [None]

    all_rows = []

    print("=== miRNA Experiment Search ===")
    print(f"Config: {config_path}")
    print(f"Sources: {', '.join(sources)}")
    print(f"miRNAs: {len(mirnas)}")
    print(f"Cell lines: {len(cell_lines)}")
    print(f"Experiment types: {as_list(config.get('experiment_types'))}")
    print(f"Perturbation types: {as_list(config.get('perturbation_types'))}")
    print(f"Species: {config.get('species')}")
    print(f"Require control: {config.get('require_control')}")

    for source, mirna, cell_line in product(sources, mirna_values, cell_line_values):
        try:
            rows = search_combination(
                source=source,
                mirna=mirna,
                cell_line=cell_line,
                config=config,
            )
            all_rows.extend(rows)
        except Exception as exc:
            print(f"WARNING: search failed for source={source}, mirna={mirna}, cell_line={cell_line}: {exc}")

    if all_rows:
        candidates = pd.DataFrame(all_rows)

        candidates = candidates.drop_duplicates(
            subset=[
                "source",
                "uid",
                "accession",
                "query_mirna",
                "query_cell_line",
            ]
        )

        candidates = apply_post_filters(candidates, config)

        candidates = candidates.sort_values(
            [
                "source",
                "query_mirna",
                "query_cell_line",
                "accession",
            ],
            ascending=True,
        )
    else:
        candidates = pd.DataFrame()

    candidate_path = outdir / "candidate_experiments.tsv"
    summary_path = outdir / "summary_statistics.tsv"
    by_mirna_path = outdir / "summary_by_mirna.tsv"
    by_accession_path = outdir / "summary_by_accession.tsv"

    candidates.to_csv(candidate_path, sep="\t", index=False)

    summary = make_summary_statistics(candidates)
    summary.to_csv(summary_path, sep="\t", index=False)

    by_mirna = make_summary_by_mirna(candidates)
    by_mirna.to_csv(by_mirna_path, sep="\t", index=False)

    by_accession = make_summary_by_accession(candidates)
    by_accession.to_csv(by_accession_path, sep="\t", index=False)

    print("\n=== Done ===")
    print(f"Candidate experiments: {len(candidates)}")
    print(f"Saved: {candidate_path}")
    print(f"Saved: {summary_path}")
    print(f"Saved: {by_mirna_path}")
    print(f"Saved: {by_accession_path}")

    if not summary.empty:
        print("\nSummary:")
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()