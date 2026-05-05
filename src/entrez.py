import time
from typing import Any

import requests


NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def entrez_esearch(
    db: str,
    query: str,
    retmax: int = 20,
    sleep_seconds: float = 0.34,
) -> list[str]:
    """
    Run NCBI Entrez esearch and return UID list.
    """
    url = f"{NCBI_EUTILS}/esearch.fcgi"

    params = {
        "db": db,
        "term": query,
        "retmode": "json",
        "retmax": retmax,
    }

    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()

    time.sleep(sleep_seconds)

    data = response.json()
    return data.get("esearchresult", {}).get("idlist", [])


def entrez_esummary(
    db: str,
    ids: list[str],
    sleep_seconds: float = 0.34,
) -> list[dict[str, Any]]:
    """
    Run NCBI Entrez esummary for a list of UIDs.
    """
    if not ids:
        return []

    url = f"{NCBI_EUTILS}/esummary.fcgi"

    params = {
        "db": db,
        "id": ",".join(ids),
        "retmode": "json",
    }

    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()

    time.sleep(sleep_seconds)

    data = response.json()
    result = data.get("result", {})

    records = []

    for uid in result.get("uids", []):
        record = result.get(uid, {})
        record["uid"] = uid
        records.append(record)

    return records


def search_entrez(
    source: str,
    query: str,
    retmax: int = 20,
    sleep_seconds: float = 0.34,
) -> list[dict[str, Any]]:
    """
    Search GEO or SRA using Entrez.

    source:
      GEO -> db='gds'
      SRA -> db='sra'
    """
    source = source.upper()

    if source == "GEO":
        db = "gds"
    elif source == "SRA":
        db = "sra"
    else:
        raise ValueError(f"Unsupported source: {source}. Use GEO or SRA.")

    ids = entrez_esearch(
        db=db,
        query=query,
        retmax=retmax,
        sleep_seconds=sleep_seconds,
    )

    records = entrez_esummary(
        db=db,
        ids=ids,
        sleep_seconds=sleep_seconds,
    )

    for record in records:
        record["source"] = source
        record["entrez_db"] = db
        record["search_query"] = query

    return records