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
`imaging/func_xcpd/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a` (functional outputs) and `XCP-D anat` is
`imaging/structure_function/hbcd/derivatives/xcp_d-0f306a2f+0ef9c88a` (anatomical outputs).

| File | Contents |
|---|---|
| `demographics_sed/sed_basic_demographics.tsv` | Infant sex, gestational age at birth, maternal age, education, income, race and ethnicity |
| `parent_visit/par_visit_data.tsv` | Recruitment site |
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
| `XCP-D anat/sub-*/ses-V02/anat/*_den-91k_{thickness,sulc}.dscalar.nii`, `*_midthickness.surf.gii` | Vertex-wise thickness and sulcal depth, midthickness surfaces |
| `XCP-D anat/sub-*/ses-V02/anat/*_seg-MyersLabonte_stat-mean_desc-{thickness,sulc}_morph.tsv` | Myers-Labonté parcel thickness and sulcal depth |
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

## Tests and formatting without HBCD data

Run `make format` to format Python with Black and R with Air. Use
`make format-python` or `make format-r` to format only one language.
Run `make test` for both the Python pytest and R testthat suites. Use
`make test-python` or `make test-r` to run one language. These targets use uv and the
formatters declared in the project's development dependencies.

```bash
uv sync --dev
uv run pytest
uv run pytest --cov --cov-report=term-missing
```

The Python suite generates entirely synthetic participant tables, NIfTI volumes,
CIFTI scalar/label files, GIFTI surfaces, connectivity matrices, and motion files
under pytest temporary directories. It does not require the controlled-access
release, atlas downloads, R, or a local `.env`. All Python functions are exercised,
including preparation and reporting entry points; R analyses are outside this
suite's scope. Fixtures are derived from the schemas used by the code, rather than
verified against an available release header export.

Numerical checks cover scanner normalization, white-matter mask resampling,
invalid imaging values, surface-area weighting, censoring, signed network metrics,
participant inclusion, units, missing-data denominators, and correlation reporting.
Integration tests reopen the generated CSV/Word files and render a synthetic
surface figure. They do not establish scientific validity on real HBCD data.

Preparation scripts can now be imported without running a pipeline. Their `main()`
functions accept `expected_count` (and, for functional data, `expected_parcels`),
so tests can validate small releases. Direct script execution retains the original
Release 2.1 count checks: 1,595 diffusion participants, 1,442 structural/functional
scans, and 276 retained functional parcels.

Coverage must remain at least 98% of Python statements. The tests call entry
points directly; five `if __name__ == "__main__"` dispatch calls are the only
currently uncovered statements.

The initial R suite in `tests/testthat/` requires `testthat`, `withr`, and `mgcv`.
If needed, install the test dependencies into your R library with:

```bash
uv run Rscript --vanilla -e 'install.packages(c("testthat", "withr", "mgcv"), repos="https://cloud.r-project.org")'
make test-r
```

R fixtures use temporary release/data directories and do not read the repository's
`.env`. Tests cover shared file I/O, participant joins, scanner-batch thresholds,
complete-case selection, formula construction, standardization, and GAM estimation.
Effect estimates and robust SEs are checked against a hand-computed example and
independent OLS/HC1 matrix algebra. Integration tests execute
`prepare/06_model_datasets.R` against invented release tables to verify visit/run
matching, age conversion, modality selection, and the five-minute retention boundary.
The covariate-preparation integration tests also run `prepare/01_covariates.R`
against an invented 2,175-person exposure cohort that satisfies the production
release-count check. They verify Release 2.1 paths, longitudinal demographic/site
selection, exposure groups, missing complications, and tissue-volume exclusions.
Python tests independently check anatomical versus functional release paths,
Pearson-based network metrics (including perfect correlations and single-sign
graphs), parcel alignment, and Figure 4 tiling and rendering.

The real CovBat test is explicitly skipped if `ComBatFamily` is not installed;
harmonization is never replaced with a mock. Run `make test-r-full` to require that
dependency and fail instead of skipping. Install ComBatFamily as described in
[its installation instructions](https://github.com/andy1764/ComBatFamily#1-installation).
Shared R utilities now load that dependency only when `harmonize()` is called.

This is an initial R suite, not full R function or pipeline coverage. It does not
yet test the full global-model pipeline,
delivery-timing bootstrap, sensitivity analyses, regional inference, or R figures.
Passing these tests does not validate the complete scientific analysis. In
particular, the ComBat test checks a balanced synthetic example; it does not establish inferential
calibration under scanner/exposure imbalance. Python coverage percentages do not
include R code.

## Extracting release headers

From the repository directory, with `HBCD_RELEASE` exported and pointing to the raw HBCD release
root (the directory containing `imaging/`, `pregnancy_exposure/`, etc.):

```bash
uv run bash scripts/extract_headers.sh > hbcd_headers.json
```

The extractor uses the existing `HBCD_RELEASE` environment variable;
`HBCD_HDP_DATA` is not needed because these are raw release files.
On a data machine without uv, `bash scripts/extract_headers.sh > hbcd_headers.json`
also works; the script requires only Bash and Python 3's standard library.
The JSON contains each distinct ordered header for the requested file types,
including hash suffixes, and ordered parcel row labels for connectivity matrices.
It records generic file patterns and missing file types, but no participant-specific
filenames or data values. Matrix row labels must match anatomical parcel naming;
unexpected labels or participant-like identifiers cause extraction to fail without
printing a partial report. Check the exit status and review the JSON before sharing.
The script does not extract the separate CIFTI atlas label mapping or data-dictionary
definitions of units, missing-value codes, and categories.

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
