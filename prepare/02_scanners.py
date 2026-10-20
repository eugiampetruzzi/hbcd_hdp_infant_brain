"""Scanner serial number and software version per V02 modality, from the BIDS scans.tsv files."""

import csv
import re
from collections import defaultdict

import pandas as pd

from config import DATA, release

FIELDS = (
    "scanner_manufacturer",
    "scanner_model",
    "scanner_software",
    "scanner_serial_number",
)
MANUFACTURERS = {
    "Siemens Healthineers": "Siemens",
    "GE MEDICAL SYSTEMS": "GE",
    "Philips Medical Systems": "Philips",
}
MODALITIES = {
    "dwi": ("dwi/", "_dwi.nii"),
    "t2w": ("anat/", "_T2w.nii"),
    "bold": ("func/", "_bold.nii"),
}


def clean(field, value):
    value = (value or "").strip()
    if field == "scanner_manufacturer":
        value = MANUFACTURERS.get(value, value)
    if field == "scanner_software":
        value = value.split("\\\\")[0]
    if field == "scanner_serial_number" and re.fullmatch(r"\d+\.0", value):
        value = value[:-2]
    return value


def main():
    scanners = defaultdict(lambda: defaultdict(set))
    for path in sorted(
        release("imaging/qc_mriqc/hbcd/rawdata").glob(
            "sub-*/ses-V02/*_ses-V02_scans.tsv"
        )
    ):
        with path.open(errors="ignore") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                name = row.get("filename", "")
                metadata = tuple(clean(field, row.get(field)) for field in FIELDS)
                for modality, (folder, suffix) in MODALITIES.items():
                    if name.startswith(folder) and suffix in name and any(metadata):
                        scanners[path.parent.parent.name][modality].add(metadata)

    rows = []
    for participant, found in sorted(scanners.items()):
        row = {"participant_id": participant}
        for modality in MODALITIES:
            values = found.get(modality, set())
            # Serial and software are reported only when every acquisition used one scanner.
            single = next(iter(values)) if len(values) == 1 else ("",) * len(FIELDS)
            row[f"{modality}_unique_scanners"] = len(values)
            row[f"{modality}_scanner_serial_number"] = single[3]
            row[f"{modality}_scanner_software"] = single[2]
        rows.append(row)

    out = DATA / "covariates/v02_scanner_metadata.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
