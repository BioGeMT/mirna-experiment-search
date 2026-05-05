import re

from src.vocab import (
    DEFAULT_EXCLUDE_TERMS,
    EXPERIMENT_TYPES,
    PERTURBATION_TYPES,
    SPECIES_TERMS,
)


def quote_term(term: str) -> str:
    term = str(term).strip()
    return f'"{term}"'


def or_group(terms: list[str]) -> str:
    clean_terms = sorted(set(str(t).strip() for t in terms if str(t).strip()))

    if not clean_terms:
        return ""

    return "(" + " OR ".join(quote_term(t) for t in clean_terms) + ")"


def and_join(groups: list[str]) -> str:
    clean_groups = [g for g in groups if g]

    if not clean_groups:
        return ""

    return " AND ".join(clean_groups)


def mirna_aliases(mirna: str) -> list[str]:
    """
    Generate safe mature-miRNA aliases.

    Examples:
      hsa-miR-22-3p -> hsa-miR-22-3p, miR-22-3p, miR-22, MIR22
      hsa-miR-129-2-3p -> hsa-miR-129-2-3p, miR-129-2-3p, miR-129-2, miR-129

    This intentionally avoids generic aliases like miR or mir.
    """
    mirna = str(mirna).strip()

    if not mirna:
        return []

    x = mirna.replace("hsa-", "")

    aliases = {mirna, x}

    # Remove mature-arm suffix: miR-22-3p -> miR-22
    base = re.sub(r"-(3p|5p)$", "", x)
    aliases.add(base)

    # Locus-style names: miR-129-2-3p -> miR-129
    match = re.match(r"^(miR-\d+[a-z]?)-\d+-(3p|5p)$", x)
    if match:
        aliases.add(match.group(1))

    expanded = set()

    for alias in aliases:
        expanded.add(alias)

        if alias.startswith("miR-"):
            expanded.add(alias.replace("miR-", "MIR"))

        if alias.startswith("miR"):
            expanded.add(alias.replace("miR", "mir", 1))

    bad_aliases = {
        "",
        "miR",
        "mir",
        "MIR",
        "miRNA",
        "microRNA",
    }

    expanded = {
        a for a in expanded
        if a not in bad_aliases and len(a) >= 5
    }

    return sorted(expanded, key=len, reverse=True)


def expand_selected_terms(
    selected: list[str] | None,
    vocabulary: dict[str, list[str]],
    label: str,
) -> list[str]:
    if not selected:
        return []

    terms: list[str] = []

    for item in selected:
        if item not in vocabulary:
            allowed = ", ".join(sorted(vocabulary.keys()))
            raise ValueError(
                f"Unknown {label}: {item}. Allowed values: {allowed}"
            )
        terms.extend(vocabulary[item])

    return terms


def expand_species(species: str | None) -> list[str]:
    if not species:
        return []

    if species in SPECIES_TERMS:
        return SPECIES_TERMS[species]

    # Allow custom species text if not built in
    return [species]


def build_query(
    mirna: str | None = None,
    cell_line: str | None = None,
    species: str | None = None,
    experiment_types: list[str] | None = None,
    perturbation_types: list[str] | None = None,
    require_control: bool = False,
    control_terms: list[str] | None = None,
    exclude_precursor_terms: bool = True,
    date_from: int | None = None,
    date_to: int | None = None,
) -> str:
    groups: list[str] = []

    if mirna:
        groups.append(or_group(mirna_aliases(mirna)))

    if cell_line:
        groups.append(or_group([cell_line]))

    species_terms = expand_species(species)
    if species_terms:
        groups.append(or_group(species_terms))

    experiment_terms = expand_selected_terms(
        experiment_types,
        EXPERIMENT_TYPES,
        label="experiment type",
    )
    if experiment_terms:
        groups.append(or_group(experiment_terms))

    perturbation_terms = expand_selected_terms(
        perturbation_types,
        PERTURBATION_TYPES,
        label="perturbation type",
    )
    if perturbation_terms:
        groups.append(or_group(perturbation_terms))

    if require_control and control_terms:
        groups.append(or_group(control_terms))

    query = and_join(groups)

    if exclude_precursor_terms:
        query += " NOT " + or_group(DEFAULT_EXCLUDE_TERMS)

    if date_from or date_to:
        start = date_from if date_from else 1900
        end = date_to if date_to else 3000
        query += f' AND ("{start}"[PDAT] : "{end}"[PDAT])'

    return query