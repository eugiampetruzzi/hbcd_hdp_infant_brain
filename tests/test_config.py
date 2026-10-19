import os

import pytest


def test_env_file_preserves_explicit_settings(modules, tmp_path, monkeypatch):
    m = modules("python/config.py")
    monkeypatch.setenv("SYNTH_EXISTING", "explicit")
    monkeypatch.delenv("SYNTH_NEW", raising=False)
    path = tmp_path / "settings.env"
    path.write_text(
        "#comment=value\nSYNTH_EXISTING=file\nSYNTH_NEW = a=b \ninvalid\n =ignored\n"
    )
    m._read_env(path)
    assert os.environ["SYNTH_EXISTING"] == "explicit"
    assert os.environ["SYNTH_NEW"] == "a=b"
    monkeypatch.delenv("SYNTH_NEW")
    m._read_env(tmp_path / "absent")


def test_release_requires_explicit_root(modules, tmp_path, monkeypatch):
    m = modules("python/config.py")
    assert m.release("imaging/file.tsv") == tmp_path / "release/imaging/file.tsv"
    monkeypatch.delenv("HBCD_RELEASE")
    with pytest.raises(SystemExit, match="HBCD_RELEASE"):
        m.release("imaging/file.tsv")


def test_data_root_required_without_env_file(modules, monkeypatch):
    # Prevent this test from reading a developer's actual .env.
    from pathlib import Path

    monkeypatch.delenv("HBCD_HDP_DATA")
    original_exists = Path.exists
    monkeypatch.setattr(
        Path, "exists", lambda p: False if p.name == ".env" else original_exists(p)
    )
    with pytest.raises(SystemExit, match="HBCD_HDP_DATA"):
        modules("python/config.py")
