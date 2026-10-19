import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _read_env(path):
    if path.exists():
        for line in path.read_text().splitlines():
            key, sep, value = line.partition("=")
            if sep and key.strip() and not key.startswith("#"):
                os.environ.setdefault(key.strip(), value.strip())


_read_env(ROOT / ".env")
if not os.environ.get("HBCD_HDP_DATA"):
    raise SystemExit("Set HBCD_HDP_DATA in .env (see .env.example)")

DATA = Path(os.environ["HBCD_HDP_DATA"])
TABLES = DATA / "derivatives" / "tables"
FIGURES = DATA / "derivatives" / "figures"
ATLAS = DATA / "atlas"

# Raw release locations, used only by prepare/. XCPD is relative to HBCD_RELEASE.
XCPD = "imaging/func_xcpd/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"


def release(rel):
    if not os.environ.get("HBCD_RELEASE"):
        raise SystemExit("Set HBCD_RELEASE in .env (see .env.example)")
    return Path(os.environ["HBCD_RELEASE"]) / rel
