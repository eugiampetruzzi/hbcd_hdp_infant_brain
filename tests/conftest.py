"""Isolated synthetic data; explicit paths override local configuration."""

import importlib.util
import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def modules(tmp_path, monkeypatch):
    monkeypatch.setenv("HBCD_HDP_DATA", str(tmp_path / "data"))
    monkeypatch.setenv("HBCD_RELEASE", str(tmp_path / "release"))
    monkeypatch.syspath_prepend(str(ROOT / "python"))
    original_exists = Path.exists
    monkeypatch.setattr(
        Path, "exists", lambda p: False if p == ROOT / ".env" else original_exists(p)
    )
    monkeypatch.delitem(sys.modules, "config", raising=False)
    import config

    def load(relative):
        spec = importlib.util.spec_from_file_location(
            Path(relative).stem, ROOT / relative
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    yield load
    sys.modules.pop("config", None)
