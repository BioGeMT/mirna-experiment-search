from pathlib import Path
from typing import Any

import yaml


def read_lines_file(path: str | None) -> list[str]:
    """
    Read a simple TXT file with one item per line.

    Empty lines and comment lines starting with # are ignored.
    """
    if path is None:
        return []

    file_path = Path(path)

    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    values: list[str] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            x = line.strip()
            if not x:
                continue
            if x.startswith("#"):
                continue
            values.append(x)

    return values


def load_yaml(path: str) -> dict[str, Any]:
    file_path = Path(path)

    if not file_path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with file_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if data is None:
        return {}

    if not isinstance(data, dict):
        raise ValueError("Config YAML must contain a dictionary at the top level.")

    return data