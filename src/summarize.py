import pandas as pd


def split_semicolon_column(df: pd.DataFrame, column: str) -> pd.Series:
    """
    Split semicolon-separated values and return exploded series.
    """
    if column not in df.columns or df.empty:
        return pd.Series(dtype="object")

    values = (
        df[column]
        .fillna("")
        .astype(str)
        .str.split(";")
        .explode()
        .str.strip()
    )

    return values[values != ""]


def count_unique_nonempty(df: pd.DataFrame, column: str) -> int:
    if column not in df.columns or df.empty:
        return 0

    return df[column].dropna().astype(str).replace("", pd.NA).dropna().nunique()


def make_summary_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create compact summary statistics for candidate experiment records.
    """
    rows = []

    def add(metric: str, count: int):
        rows.append({"metric": metric, "count": int(count)})

    add("total_candidate_records", len(df))

    if df.empty:
        return pd.DataFrame(rows)

    add("unique_accessions", count_unique_nonempty(df, "accession"))
    add("unique_query_mirnas", count_unique_nonempty(df, "query_mirna"))
    add("unique_query_cell_lines", count_unique_nonempty(df, "query_cell_line"))

    if "source" in df.columns:
        for source, n in df["source"].fillna("unknown").value_counts().items():
            add(f"records_source_{source}", n)

    if "has_control_keyword" in df.columns:
        add("records_with_control_keyword", df["has_control_keyword"].sum())
        add("records_without_control_keyword", (~df["has_control_keyword"]).sum())

    if "has_perturbation_keyword" in df.columns:
        add("records_with_perturbation_keyword", df["has_perturbation_keyword"].sum())

    if "has_overexpression_keyword" in df.columns:
        add("records_with_overexpression_keyword", df["has_overexpression_keyword"].sum())

    if "has_knockdown_keyword" in df.columns:
        add("records_with_knockdown_keyword", df["has_knockdown_keyword"].sum())

    if "has_knockout_keyword" in df.columns:
        add("records_with_knockout_keyword", df["has_knockout_keyword"].sum())

    experiment_values = split_semicolon_column(df, "matched_experiment_types")
    for experiment_type, n in experiment_values.value_counts().items():
        add(f"records_experiment_type_{experiment_type}", n)

    perturbation_values = split_semicolon_column(df, "matched_perturbation_types")
    for perturbation_type, n in perturbation_values.value_counts().items():
        add(f"records_perturbation_type_{perturbation_type}", n)

    cell_line_values = split_semicolon_column(df, "inferred_cell_lines")
    for cell_line, n in cell_line_values.value_counts().items():
        add(f"records_inferred_cell_line_{cell_line}", n)

    return pd.DataFrame(rows)


def make_summary_by_mirna(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "query_mirna" not in df.columns:
        return pd.DataFrame()

    return (
        df.groupby("query_mirna", dropna=False)
        .agg(
            n_records=("query_mirna", "size"),
            n_unique_accessions=("accession", "nunique"),
            n_geo_records=("source", lambda x: (x == "GEO").sum()),
            n_sra_records=("source", lambda x: (x == "SRA").sum()),
            n_with_control=("has_control_keyword", "sum"),
            inferred_cell_lines=(
                "inferred_cell_lines",
                lambda x: ";".join(sorted(set(v for s in x.dropna().astype(str) for v in s.split(";") if v))),
            ),
            matched_experiment_types=(
                "matched_experiment_types",
                lambda x: ";".join(sorted(set(v for s in x.dropna().astype(str) for v in s.split(";") if v))),
            ),
            matched_perturbation_types=(
                "matched_perturbation_types",
                lambda x: ";".join(sorted(set(v for s in x.dropna().astype(str) for v in s.split(";") if v))),
            ),
        )
        .reset_index()
        .sort_values("n_records", ascending=False)
    )


def make_summary_by_accession(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "accession" not in df.columns:
        return pd.DataFrame()

    return (
        df.groupby(["source", "accession"], dropna=False)
        .agg(
            n_records=("accession", "size"),
            query_mirnas=(
                "query_mirna",
                lambda x: ";".join(sorted(set(v for v in x.dropna().astype(str) if v))),
            ),
            query_cell_lines=(
                "query_cell_line",
                lambda x: ";".join(sorted(set(v for v in x.dropna().astype(str) if v))),
            ),
            has_control_keyword=("has_control_keyword", "max"),
            inferred_cell_lines=(
                "inferred_cell_lines",
                lambda x: ";".join(sorted(set(v for s in x.dropna().astype(str) for v in s.split(";") if v))),
            ),
            matched_experiment_types=(
                "matched_experiment_types",
                lambda x: ";".join(sorted(set(v for s in x.dropna().astype(str) for v in s.split(";") if v))),
            ),
            matched_perturbation_types=(
                "matched_perturbation_types",
                lambda x: ";".join(sorted(set(v for s in x.dropna().astype(str) for v in s.split(";") if v))),
            ),
            title=("title", "first"),
            url=("url", "first"),
        )
        .reset_index()
        .sort_values("n_records", ascending=False)
    )