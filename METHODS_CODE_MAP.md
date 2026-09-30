# Methods to code map

| Manuscript component | Reproducible code |
|---|---|
| HDP and control definitions | `derive/01_build_exposure.py` |
| Diffusion MRI outcomes | `derive/02_diffusion_global.py`, `derive/03_diffusion_bundles.py` |
| Structural MRI outcomes | `derive/04_structural_morphometry.py` |
| Functional MRI outcomes | `derive/05_functional_networks.py` |
| Imaging and covariate merges | `derive/06_build_analysis_datasets.py` |
| Primary global HDP analyses | `analysis/01_primary_global.R` |
| Gestational-age pathway | `analysis/02_gestational_mediation.R` |
| Other-complication sensitivity | `analysis/03_no_complications.R` |
| HDP clinical presentations | `analysis/04_hdp_presentations.R` |
| Scanner and fMRI thresholds | `analysis/05_mri_sensitivity.R` |
| Maternal covariates | `analysis/06_covariate_sensitivity.R` |
| Spatial diffusion analyses | `analysis/07_spatial_diffusion.R` |
| Spatial structural analyses | `analysis/08_spatial_structural.R` |
| Spatial functional analyses | `analysis/09_spatial_functional.R` |
| Tables 1 and 2 | `analysis/10_descriptive_tables.R`, `tables/01_build_tables.py` |
| Technical validation | `validation/01_validate_results.py` |

The current primary control definition is applied centrally in `analysis/common.R`: controls must have `primary_hdp == 0` and `none_of_these == 1`; all eligible HDP participants are retained regardless of other complications. The additional sensitivity analysis requires `other_complication == 0` among HDP participants.
