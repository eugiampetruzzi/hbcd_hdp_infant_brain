#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

from hbcd_hdp.config import get_settings
from hbcd_hdp.io import sha256, write_tsv


def main() -> None:
    if sys.version_info[:2] != (3, 11):
        raise SystemExit(f"Python 3.11 required; found {sys.version.split()[0]}")
    settings = get_settings()
    manifest_path = Path(__file__).parents[1] / "config/input_manifest.tsv"
    with manifest_path.open() as stream:
        manifest = list(csv.DictReader(stream, delimiter="\t"))
    rows = []
    missing = []
    for item in manifest:
        path = settings.root / item["path_relative_to_HBCD_ROOT"]
        exists = path.is_file()
        rows.append(
            {
                **item,
                "absolute_path": str(path),
                "exists": int(exists),
                "bytes": path.stat().st_size if exists else "",
                "sha256": sha256(path) if exists else "",
            }
        )
        if item["required"] == "yes" and not exists:
            missing.append(path)
    report = settings.output / "provenance/input_manifest_resolved.tsv"
    write_tsv(__import__("pandas").DataFrame(rows), report)
    if missing:
        raise SystemExit("Missing required inputs:\n" + "\n".join(map(str, missing)))
    print(f"Validated {len(rows)} inputs; wrote {report}")


if __name__ == "__main__":
    main()
