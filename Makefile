-include .env
export

PYTHON ?= uv run python

.PHONY: all prepare models timing sensitivity regional figures tables clean format format-python format-r test test-python test-r test-r-full

all: models timing sensitivity regional figures tables

format: format-python format-r

format-python:
	uv run black python prepare tests

format-r:
	uv run air format R prepare tests/testthat tests/testthat.R

test: test-python test-r

test-python:
	uv run pytest

test-r:
	uv run Rscript --vanilla tests/testthat.R

test-r-full:
	uv run Rscript --vanilla tests/testthat.R --require-combat

# Builds the derived inputs from the raw release (HBCD_RELEASE); run on Carina.
prepare:
	Rscript prepare/01_covariates.R
	PYTHONPATH=python $(PYTHON) prepare/02_scanners.py
	PYTHONPATH=python $(PYTHON) prepare/03_diffusion.py
	PYTHONPATH=python $(PYTHON) prepare/04_structural.py
	PYTHONPATH=python $(PYTHON) prepare/05_functional.py
	Rscript prepare/06_model_datasets.R

models:
	Rscript R/01_global_models.R

timing: models
	Rscript R/02_delivery_timing.R

sensitivity: models
	Rscript R/03_sensitivity.R

regional: models
	Rscript R/04_regional.R

figures: timing regional
	Rscript R/05_figures.R
	$(PYTHON) python/figure4_map.py

tables: sensitivity
	$(PYTHON) python/tables.py

clean:
	@test -n "$(HBCD_HDP_DATA)" || (echo "Set HBCD_HDP_DATA in .env" && exit 1)
	rm -rf "$(HBCD_HDP_DATA)/derivatives"
