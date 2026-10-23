"""Synthetic records using the supplied release's actual column names."""

import json
from functools import lru_cache
from pathlib import Path

import pandas as pd


@lru_cache
def release_schemas():
    return json.loads(
        (Path(__file__).resolve().parents[1] / "hbcd_headers.json").read_text()
    )["files"]


def release_columns(kind, variant=0):
    return list(release_schemas()[kind]["schemas"][variant]["columns"])


def release_frame(kind, rows, variant=0):
    columns = release_columns(kind, variant)
    frame = pd.DataFrame(rows)
    assert set(frame.columns) <= set(columns)
    return frame.reindex(columns=columns)
