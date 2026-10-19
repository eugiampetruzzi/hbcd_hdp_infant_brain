"""Global cerebral white-matter FA, MD and MKT, and bundle medians, from QSIRecon DIPYDKI."""

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import affine_transform

from config import DATA, release

QSIRECON = release("imaging/dwi_qsirecon")
QSIPREP = release("imaging/dwi_qsiprep/hbcd/derivatives/qsiprep")
MAPS = {"mkt": "model-dki_param-mkt", "fa": "model-tensor_param-fa", "md": "model-dki_param-md"}
WM_LABELS = (2, 41)


# The release is split across several qsirecon-DIPYDKI folders; where a participant has
# copies in more than one, the last in path order is used.
def by_participant(suffix):
    paths = sorted(QSIRECON.glob(f"**/qsirecon-DIPYDKI/sub-*/ses-V02/dwi/*_ses-V02_space-ACPC_{suffix}"))
    return {p.name.split("_")[0]: p for p in paths}


def white_matter_medians(participant, maps):
    grid = nib.load(maps["mkt"])
    name = f"{participant}_ses-V02_space-ACPC_desc-aseg_dseg.nii.gz"
    dseg = nib.load(QSIPREP / participant / "ses-V02/anat" / name)
    to_dseg = np.linalg.inv(dseg.affine) @ grid.affine
    labels = affine_transform(np.asarray(dseg.dataobj), to_dseg[:3, :3], offset=to_dseg[:3, 3],
                              output_shape=grid.shape, order=0, mode="constant", cval=0, prefilter=False)
    # Infants missing either hemisphere's label are excluded (incomplete segmentation).
    if not all((labels == label).any() for label in WM_LABELS):
        return {f"global_median_{m}": np.nan for m in MAPS}
    wm = np.isin(labels, WM_LABELS)
    row = {}
    for metric, path in maps.items():
        values = nib.load(path).get_fdata()[wm]
        values = values[np.isfinite(values) & (values != 0)]
        row[f"global_median_{metric}"] = np.median(values) if values.size else np.nan
    return row


maps = {m: by_participant(f"{entities}_dwimap.nii.gz") for m, entities in MAPS.items()}
participants = sorted(set.intersection(*(set(v) for v in maps.values())))
assert len(participants) == 1595, len(participants)
rows = [{"participant_id": p, "session_id": "ses-V02",
         **white_matter_medians(p, {m: maps[m][p] for m in MAPS})} for p in participants]
out = DATA / "outcomes/diffusion"
out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out / "v02_global_cerebral_wm_fa_md.tsv", sep="\t", index=False)

columns = ["subject_id", "bundle", "variable_name", "masked_median"]
stats = pd.concat(pd.read_csv(p, sep="\t", usecols=columns)
                  for p in by_participant("bundles-DSIStudio_scalarstats.tsv").values())
stats["masked_median"] = stats.masked_median.replace([np.inf, -np.inf], np.nan)
for metric in MAPS:
    wide = stats[stats.variable_name == f"dki_{metric}"].pivot(index="subject_id", columns="bundle",
                                                                 values="masked_median")
    wide = wide.add_prefix(f"{metric}_").rename_axis(index="participant_id", columns=None).reset_index()
    wide.insert(1, "session_id", "ses-V02")
    wide.to_csv(out / f"v02_bundle_{metric}_masked_median.tsv", sep="\t", index=False)
