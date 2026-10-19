# hbcd_hdp_infant_brain

Analysis code for *A cortical signature of hypertensive pregnancy in the infant brain*. Using
multimodal MRI from the HEALthy Brain and Child Development (HBCD) Study (Release 2.1), the
pipeline harmonizes global and regional imaging measures across scanners with CovBat-GAM, compares
infants born after hypertensive and normotensive pregnancies at matched postmenstrual age, and
examines how delivery timing relates to these differences.

## Data access

HBCD data are access-controlled and are not included here. An approved HBCD Data Use Certification
is required (https://hbcdstudy.org). The preparation stage reads the Release 2.1 files below from the
release root (`HBCD_RELEASE`) and writes the derived inputs to a data directory outside this
repository (`HBCD_HDP_DATA`). The release is stored on Carina, so `make prepare` must be run there.
All analysis outputs are written to `derivatives/` inside the data directory.

### Raw HBCD Release 2.1 files

Paths are relative to the release root. `XCP-D` is
`imaging/func_xcpd/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a`.

| File | Contents |
|---|---|
| `demographics/sed_basic_demographics.tsv` | Infant sex, gestational age at birth, maternal age, education, income, race and ethnicity |
| `demographics/par_visit_data.tsv` | Recruitment site |
| `pregnancy_exposure/pex_bm_healthv2_preg.tsv` | Visit 2 pregnancy health (complications, aspirin, labor, delivery, medications) |
| `pregnancy_exposure/pex_bm_health_preg__chroncond.tsv` | Visit 1 chronic conditions (hypertension, diabetes) |
| `pregnancy_exposure/pex_bm_healthv2_inf.tsv` | Visit 2 infant health (birth weight, NICU) |
| `pregnancy_exposure/pex_bm_health_preg__{healthhx,meds}.tsv`, `pex_ch_tlfb.tsv` | Visit 1 health history and medications; Timeline Followback |
| `imaging/qc_mriqc/hbcd/rawdata/sub-*/ses-V02/*_scans.tsv` | Scanner serial number and software version per acquisition |
| `imaging/qc_brainswipes/img_brainswipes_xcpd_hash-0f306a2f+0ef9c88a_T2w.tsv` | BrainSwipes score and postnatal age at the structural scan |
| `imaging/struct_bibsnet/img_bibsnet_space-T2w_desc-aseg_volumes.tsv` | BIBSNet T2w segmentation volumes |
| `imaging/dwi_qsiprep/img_qsiprep_space-ACPC_desc-image_qc.tsv` | Neighboring-DWI correlation and postnatal age at the diffusion scan |
| `imaging/dwi_qsiprep/hbcd/derivatives/qsiprep/sub-*/ses-V02/anat/*_space-ACPC_desc-aseg_dseg.nii.gz` | SynthSeg segmentation (cerebral white matter) |
| `imaging/dwi_qsirecon/**/qsirecon-DIPYDKI/sub-*/ses-V02/dwi/*_dwimap.nii.gz` | FA (tensor fit), MD and MKT (kurtosis fit) maps |
| `imaging/dwi_qsirecon/**/qsirecon-DIPYDKI/sub-*/ses-V02/dwi/*_bundles-DSIStudio_scalarstats.tsv` | Bundle FA, MD and MKT medians |
| `XCP-D/sub-*/ses-V02/anat/*_den-91k_{thickness,sulc}.dscalar.nii`, `*_midthickness.surf.gii` | Vertex-wise thickness and sulcal depth, midthickness surfaces |
| `XCP-D/sub-*/ses-V02/anat/*_seg-MyersLabonte_stat-mean_desc-{thickness,sulc}_morph.tsv` | Myers-Labonté parcel thickness and sulcal depth |
| `XCP-D/sub-*/ses-V02/func/*task-rest_*_seg-MyersLabonte_stat-pearsoncorrelation_relmat.tsv` | Parcel connectivity of the concatenated resting-state runs |
| `XCP-D/sub-*/ses-V02/func/*task-rest_{motion,outliers}.tsv` | Framewise displacement and censoring mask |

### Derived inputs

The preparation stage (`prepare/`) builds these files from the raw release files, except the atlas
files, which are downloaded from the sources listed. The analysis stage reads them.

| File | Contents |
|---|---|
| `covariates/hdp_missingness_person_level.csv` | Demographics, gestational age, recruitment site |
| `covariates/v02_pregnancy_complications.tsv` | HDP exposure group and Visit 2 pregnancy-complication items |
| `covariates/v02_scanner_metadata.tsv` | Scanner serial number and software version per modality |
| `covariates/v02_dwi_primary_covariates.tsv` | Diffusion sample, exposure, covariates, neighboring-DWI correlation |
| `covariates/v02_functional_qc.tsv` | Functional motion summaries |
| `covariates/v02_bibsnet_brain_volume.tsv` | Total brain volume from BIBSNet |
| `analysis_datasets/v02_global_cortical_morphometry_model_data.tsv` | Structural sample with global thickness and sulcal depth |
| `analysis_datasets/v02_global_functional_hdp_data.tsv` | Functional sample with modularity and participation coefficient |
| `outcomes/diffusion/v02_global_cerebral_wm_fa_md.tsv` | Global white-matter FA, MD, MKT |
| `outcomes/diffusion/v02_bundle_{fa,md,mkt}_masked_median.tsv` | Bundle medians from QSIRecon |
| `outcomes/functional/v02_functional_network_metrics.tsv` | Functional network summaries |
| `outcomes/structural/v02_myerslabonte_{thickness,sulc}.tsv` | Myers-Labonté parcel thickness and sulcal depth |
| `phenotype/pex_bm_healthv2_inf.csv` | Visit 2 infant health questionnaire (birth weight, NICU) |
| `phenotype/pex_bm_healthv2_preg.csv` | Visit 2 pregnancy health (gestational diabetes, aspirin, labor, delivery, medications) |
| `phenotype/pex_bm_health_preg__healthhx.csv` | Visit 1 health history (pre-pregnancy BMI, live births, aspirin) |
| `phenotype/pex_bm_health_preg__chroncond.csv` | Visit 1 chronic conditions (diabetes) |
| `phenotype/pex_bm_health_preg__meds.csv` | Visit 1 prescribed medications (RxNorm codes and names) |
| `phenotype/pex_ch_tlfb.csv` | Timeline Followback prenatal substance use (Visits 1 and 2) |
| `atlas/parcels.dlabel.nii` | Myers-Labonté neonatal parcellation (github.com/myersm0/myers-labonte_parcellation) |
| `atlas/fsLR_{L,R}_inflated.surf.gii` | fsLR 32k inflated surfaces from TemplateFlow (`tpl-fsLR_den-32k_hemi-*_inflated.surf.gii`) |

## Setup

Requires R (4.5) with `mgcv`, `ComBatFamily` (github.com/andy1764/ComBatFamily), `sandwich`,
`ggplot2` and `patchwork`, and Python 3.11 managed with [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/eugiampetruzzi/hbcd_hdp_infant_brain
cd hbcd_hdp_infant_brain
uv sync
cp .env.example .env    # then set HBCD_HDP_DATA (and HBCD_RELEASE on Carina)
make prepare            # on Carina only
make all
```

## Pipeline

| Step | Script | Output (`derivatives/`) | Manuscript |
|---|---|---|---|
| Preparation (`make prepare`, Carina) | `prepare/01_covariates.R`, `02_scanners.py`, `03_diffusion.py`, `04_structural.py`, `05_functional.py`, `06_model_datasets.R` | Derived inputs above (data directory) | Methods |
| Harmonization and global models | `R/01_global_models.R` | `tables/global_models.csv`, `analytic_sample.csv` | Table 3 |
| Delivery timing | `R/02_delivery_timing.R` | `tables/delivery_timing_*.csv`, `term_contrast.csv`, `thickness_by_ga_band.csv`, `ga_portion_bootstrap.csv` | Figure 3, Results |
| Sensitivity analyses | `R/03_sensitivity.R` | `tables/sensitivity.csv`, `leave_one_scanner_out.csv`, `scanner_interaction.csv`, `labor_onset.csv` | Table 4 |
| Regional analyses | `R/04_regional.R` | `tables/regional.csv`, `parcel_networks.csv`, `parcel_artifact_check.csv` | Figure 4, Results |
| Figures | `R/05_figures.R`, `python/figure4_map.py` | `figures/figure{2,3,4}_*` | Figures 2-4 |
| Tables | `python/tables.py` | `tables/tables.docx` | Tables 1-4 |

`R/common.R` holds the shared steps: loading the modality samples, defining scanner batches
(serial number x software version), CovBat-GAM harmonization, and model fitting with
heteroskedasticity-robust standard errors. `analytic_sample.csv` contains participant-level data
and must not be shared.

Code was developed by E.G. with the assistance of AI coding tools (Codex, OpenAI; Claude Code,
Anthropic) and reviewed by a second author.
