#!/usr/bin/env python3

import argparse
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Search GEO/SRA for miRNA-related experiments using configurable filters."
    )

    parser.add_argument(
        "--config",
        default="config/search.yaml",
        help="Path to YAML config file.",
    )

    parser.add_argument(
        "--mirnas-file",
        default=None,
        help="Optional TXT file with one miRNA per line. Overrides config.",
    )

    parser.add_argument(
        "--cell-lines-file",
        default=None,
        help="Optional TXT file with one cell line per line. Overrides config.",
    )

    parser.add_argument(
        "--sources",
        nargs="+",
        default=None,
        choices=["GEO", "SRA"],
        help="Sources to search. Overrides config.",
    )

    parser.add_argument(
        "--species",
        default=None,
        help='Species filter, e.g. "Homo sapiens". Overrides config.',
    )

    parser.add_argument(
        "--experiment-types",
        nargs="+",
        default=None,
        help="Experiment type filters, e.g. RNA-seq CLIP CLASH qPCR. Overrides config.",
    )

    parser.add_argument(
        "--perturbation-types",
        nargs="+",
        default=None,
        help="Perturbation filters, e.g. overexpression knockdown knockout. Overrides config.",
    )

    parser.add_argument(
        "--require-control",
        action="store_true",
        help="Require control-like keyword in metadata.",
    )

    parser.add_argument(
        "--no-require-control",
        action="store_true",
        help="Disable require_control even if enabled in config.",
    )

    parser.add_argument(
        "--date-from",
        type=int,
        default=None,
        help="Publication date lower bound year. Overrides config.",
    )

    parser.add_argument(
        "--date-to",
        type=int,
        default=None,
        help="Publication date upper bound year. Overrides config.",
    )

    parser.add_argument(
        "--retmax",
        type=int,
        default=None,
        help="Max Entrez records per query/source. Overrides config.",
    )

    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=None,
        help="Delay between Entrez requests. Overrides config.",
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Overrides config.",
    )

    parser.add_argument(
        "--include-precursor-terms",
        action="store_true",
        help="Do not exclude precursor/primary miRNA terms.",
    )

    return parser.parse_args()


def apply_cli_overrides(config, args):
    config = dict(config)

    if args.mirnas_file is not None:
        config["mirnas_file"] = args.mirnas_file

    if args.cell_lines_file is not None:
        config["cell_lines_file"] = args.cell_lines_file

    if args.sources is not None:
        config["sources"] = args.sources

    if args.species is not None:
        config["species"] = args.species

    if args.experiment_types is not None:
        config["experiment_types"] = args.experiment_types

    if args.perturbation_types is not None:
        config["perturbation_types"] = args.perturbation_types

    if args.require_control:
        config["require_control"] = True

    if args.no_require_control:
        config["require_control"] = False

    if args.date_from is not None:
        config["date_from"] = args.date_from

    if args.date_to is not None:
        config["date_to"] = args.date_to

    if args.retmax is not None:
        config["retmax"] = args.retmax

    if args.sleep_seconds is not None:
        config["sleep_seconds"] = args.sleep_seconds

    if args.outdir is not None:
        config["outdir"] = args.outdir

    if args.include_precursor_terms:
        config["exclude_precursor_terms"] = False

    return config


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
    args = parse_args()

    config = load_yaml(args.config)
    config = apply_cli_overrides(config, args)

    outdir = Path(config.get("outdir", "results"))
    outdir.mkdir(parents=True, exist_ok=True)

    mirnas = read_lines_file(config.get("mirnas_file"))
    cell_lines = read_lines_file(config.get("cell_lines_file"))

    sources = as_list(config.get("sources")) or ["GEO", "SRA"]

    mirna_values = mirnas if mirnas else [None]
    cell_line_values = cell_lines if cell_lines else [None]

    all_rows = []

    print("=== miRNA Experiment Search ===")
    print(f"Config: {args.config}")
    print(f"Sources: {', '.join(sources)}")
    print(f"miRNAs: {len(mirnas)}")
    print(f"Cell lines: {len(cell_lines)}")
    print(f"Experiment types: {as_list(config.get('experiment_types'))}")
    print(f"Perturbation types: {as_list(config.get('perturbation_types'))}")
    print(f"Species: {config.get('species')}")
    print(f"Require control: {config.get('require_control')}")
    print(f"Date from: {config.get('date_from')}")
    print(f"Date to: {config.get('date_to')}")
    print(f"Outdir: {outdir}")

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
            print(
                f"WARNING: search failed for "
                f"source={source}, mirna={mirna}, cell_line={cell_line}: {exc}"
            )

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