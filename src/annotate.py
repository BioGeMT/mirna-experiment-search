import re
from typing import Any

from src.vocab import (
    COMMON_CELL_LINES,
    CONTROL_TERMS,
    EXPERIMENT_TYPES,
    PERTURBATION_TYPES,
)


def row_text(row: dict[str, Any]) -> str:
    """
    Combine likely metadata fields into one lowercase searchable string.
    """
    fields = [
        "title",
        "summary",
        "organism",
        "taxon",
        "library_strategy",
        "library_source",
        "library_selection",
        "platform",
        "samples",
        "runs",
        "experiment",
        "study",
        "bioproject",
        "biosample",
        "search_query",
    ]

    return " ".join(str(row.get(f, "")) for f in fields).lower()


def contains_term(text: str, term: str) -> bool:
    """
    Case-insensitive loose phrase matching.
    """
    term = str(term).lower().strip()
    if not term:
        return False

    pattern = r"(?<![a-z0-9])" + re.escape(term) + r"(?![a-z0-9])"
    return re.search(pattern, text) is not None


def matched_categories(
    text: str,
    vocab: dict[str, list[str]],
) -> list[str]:
    matched = []

    for category, terms in vocab.items():
        if any(contains_term(text, term) for term in terms):
            matched.append(category)

    return sorted(set(matched))


def infer_cell_lines(text: str) -> list[str]:
    found = []

    for cell_line, terms in COMMON_CELL_LINES.items():
        if any(contains_term(text, term.lower()) for term in terms):
            found.append(cell_line)

    return sorted(set(found))


def infer_has_control(text: str) -> bool:
    return any(contains_term(text, term) for term in CONTROL_TERMS)


def make_url(source: str, accession: str | None, uid: str | None) -> str:
    source = str(source).upper()

    if accession and accession != "nan":
        if source == "GEO":
            return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
        if source == "SRA":
            return f"https://www.ncbi.nlm.nih.gov/sra/?term={accession}"

    if uid and uid != "nan":
        if source == "GEO":
            return f"https://www.ncbi.nlm.nih.gov/gds/?term={uid}"
        if source == "SRA":
            return f"https://www.ncbi.nlm.nih.gov/sra/?term={uid}"

    return ""


def extract_accession(record: dict[str, Any]) -> str:
    """
    Entrez GEO and SRA summaries expose accessions differently.
    """
    for key in ["accession", "gse", "study", "experiment"]:
        value = record.get(key)
        if value:
            return str(value)

    return ""


def normalize_record(
    record: dict[str, Any],
    query_mirna: str | None = None,
    query_cell_line: str | None = None,
) -> dict[str, Any]:
    """
    Convert one raw Entrez record into a clean candidate row.
    """
    source = str(record.get("source", ""))
    uid = str(record.get("uid", ""))

    accession = extract_accession(record)

    title = record.get("title", "")
    summary = record.get("summary", "")

    organism = (
        record.get("organism")
        or record.get("taxon")
        or ""
    )

    text = row_text(record)

    experiment_matches = matched_categories(text, EXPERIMENT_TYPES)
    perturbation_matches = matched_categories(text, PERTURBATION_TYPES)
    cell_lines = infer_cell_lines(text)
    has_control = infer_has_control(text)

    return {
        "query_mirna": query_mirna or "",
        "query_cell_line": query_cell_line or "",
        "source": source,
        "uid": uid,
        "accession": accession,
        "title": title,
        "summary": summary,
        "species_or_organism": organism,
        "matched_experiment_types": ";".join(experiment_matches),
        "matched_perturbation_types": ";".join(perturbation_matches),
        "has_perturbation_keyword": bool(perturbation_matches),
        "has_overexpression_keyword": "overexpression" in perturbation_matches,
        "has_knockdown_keyword": "knockdown" in perturbation_matches,
        "has_knockout_keyword": "knockout" in perturbation_matches,
        "has_control_keyword": has_control,
        "inferred_cell_lines": ";".join(cell_lines) if cell_lines else "",
        "url": make_url(source, accession, uid),
        "search_query": record.get("search_query", ""),
    }


def normalize_records(
    records: list[dict[str, Any]],
    query_mirna: str | None = None,
    query_cell_line: str | None = None,
) -> list[dict[str, Any]]:
    return [
        normalize_record(
            record,
            query_mirna=query_mirna,
            query_cell_line=query_cell_line,
        )
        for record in records
    ]