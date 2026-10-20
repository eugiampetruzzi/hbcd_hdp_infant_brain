#!/usr/bin/env bash
# Run from the repository: datadir=/path/to/release bash scripts/extract_headers.sh > headers.json
# Requires Python 3 (standard library only). Reads source files without modifying them.
set -euo pipefail

if [[ -z "${datadir:-}" || ! -d "$datadir" ]]; then
    echo 'Set datadir to the HBCD release root directory.' >&2
    exit 1
fi
export datadir

python3 - <<'PY'
import csv
import json
import os
from pathlib import Path
import re
import sys


ANAT = "imaging/structure_function/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"
FUNC = "imaging/func_xcpd/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a"
PATTERNS = {
    "parcel_thickness": f"{ANAT}/sub-*/ses-V02/anat/*_seg-MyersLabonte_stat-mean_desc-thickness_morph.tsv",
    "parcel_sulc": f"{ANAT}/sub-*/ses-V02/anat/*_seg-MyersLabonte_stat-mean_desc-sulc_morph.tsv",
    "connectivity": f"{FUNC}/sub-*/ses-V02/func/*task-rest_space-fsLR_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv",
    "motion": f"{FUNC}/sub-*/ses-V02/func/*task-rest_motion.tsv",
    "outliers": f"{FUNC}/sub-*/ses-V02/func/*task-rest_outliers.tsv",
    "demographics": "demographics_sed/sed_basic_demographics.tsv",
    "visits": "parent_visit/par_visit_data.tsv",
    "pregnancy": "pregnancy_exposure/pex_bm_healthv2_preg.tsv",
    "chronic_conditions": "pregnancy_exposure/pex_bm_health_preg__chroncond.tsv",
    "infant_health": "pregnancy_exposure/pex_bm_healthv2_inf.tsv",
    "health_history": "pregnancy_exposure/pex_bm_health_preg__healthhx.tsv",
    "medications": "pregnancy_exposure/pex_bm_health_preg__meds.tsv",
    "substance_use": "pregnancy_exposure/pex_ch_tlfb.tsv",
    "brain_volumes": "imaging/struct_bibsnet/img_bibsnet_space-T2w_desc-aseg_volumes.tsv",
    "diffusion_qc": "imaging/dwi_qsiprep/img_qsiprep_space-ACPC_desc-image_qc.tsv",
    "structural_qc": "imaging/qc_brainswipes/img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w.tsv",
    "scanner_metadata": "imaging/qc_mriqc/hbcd/rawdata/sub-*/ses-V02/*_ses-V02_scans.tsv",
    "bundle_statistics": "imaging/dwi_qsirecon/**/qsirecon-DIPYDKI/sub-*/ses-V02/dwi/*_bundles-DSIStudio_scalarstats.tsv",
}


def parcel_label(value):
    # Only anatomical labels are exported from matrix rows, never subject labels.
    base = re.sub(r"_hash-.*$", "", value)
    return base == "None" or re.fullmatch(r"[LR]_[A-Za-z0-9_+.-]+", base)


def read_schema(path, matrix=False):
    with path.open(encoding="utf-8-sig", newline="") as source:
        rows = csv.reader(source, delimiter="\t", strict=True)
        columns = next(rows, None)
        if not columns or not any(columns):
            raise ValueError("empty header")
        if any(re.search(r"\bsub-[A-Za-z0-9]+", name) for name in columns):
            raise ValueError("participant-like identifier in header")
        schema = {"columns": columns}
        if matrix:
            if not all(parcel_label(label) for label in columns[1:]):
                raise ValueError("unrecognized matrix column labels")
            labels = []
            for row in rows:
                if not row:
                    continue
                if not parcel_label(row[0]):
                    raise ValueError("unrecognized matrix row label")
                labels.append(row[0])
            schema["row_labels"] = labels
        return schema


def collect(root, pattern, matrix):
    variants = {}
    for path in sorted(root.glob(pattern)):
        if not path.is_file():
            continue
        schema = read_schema(path, matrix)
        key = json.dumps(schema, ensure_ascii=False, sort_keys=True)
        variants[key] = schema
    return {
        "pattern": pattern,
        "status": "found" if variants else "missing",
        "schemas": [variants[key] for key in sorted(variants)],
    }


def main():
    root = Path(os.environ["datadir"])
    report = {"format_version": 1, "files": {}}
    for name, pattern in PATTERNS.items():
        try:
            report["files"][name] = collect(root, pattern, name == "connectivity")
        except (OSError, UnicodeError, csv.Error, ValueError):
            # Do not expose source paths, raw contents, or exception details.
            print(f"Cannot extract {name}: unreadable/empty TSV or unexpected identifier/parcel labels.", file=sys.stderr)
            return 1
    # Buffer until all files have been checked, so errors never produce partial JSON.
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


sys.exit(main())
PY
