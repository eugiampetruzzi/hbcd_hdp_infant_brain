SHELL := /bin/bash
PYTHON := uv run python
export PYTHONPATH := src

.PHONY: check-environment imaging datasets primary mediation secondary spatial descriptives analysis tables validate test all

check-environment:
	$(PYTHON) derive/00_validate_inputs.py

imaging: check-environment
	$(PYTHON) derive/01_build_exposure.py
	$(PYTHON) derive/02_diffusion_global.py
	$(PYTHON) derive/03_diffusion_bundles.py
	$(PYTHON) derive/04_structural_morphometry.py
	$(PYTHON) derive/05_functional_networks.py
	$(PYTHON) derive/06_build_analysis_datasets.py

datasets:
	$(PYTHON) derive/06_build_analysis_datasets.py

primary:
	Rscript analysis/01_primary_global.R

mediation:
	Rscript analysis/02_gestational_mediation.R

secondary:
	Rscript analysis/03_no_complications.R
	Rscript analysis/04_hdp_presentations.R
	Rscript analysis/05_mri_sensitivity.R
	Rscript analysis/06_covariate_sensitivity.R

spatial:
	Rscript analysis/07_spatial_diffusion.R
	Rscript analysis/08_spatial_structural.R
	Rscript analysis/09_spatial_functional.R

descriptives:
	Rscript analysis/10_descriptive_tables.R

analysis: primary mediation secondary spatial descriptives

tables:
	$(PYTHON) tables/01_build_tables.py

validate:
	$(PYTHON) validation/01_validate_results.py

test:
	uv run pytest
	Rscript tests/test_hc3.R

all: imaging analysis tables validate
