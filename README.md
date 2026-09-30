# HBCD HDP infant brain development

Code for the HBCD Release 2.1 HDP and infant MRI analyses.

- Imaging outcomes: Python 3.11
- Statistical analyses: R 4.5.1
- Carina project root: `/projects/iang/eanlp/eu_hbcd-release-2.1`

## Setup on Carina

```bash
cp .env.example .env
module load python/3.11
module load R/4.5.1
uv sync --frozen
```

## Run

```bash
make check-environment
make analysis
make tables
make validate
```

`make imaging` derives the MRI outcomes. `make all` runs the whole pipeline.

## Folders

- `derive/`: cohort construction, MRI outcome derivation, and analysis datasets
- `analysis/`: primary, mediation, secondary, spatial, and descriptive analyses
- `tables/`: formatted table exports
- `validation/`: sample-size and result checks
- `config/`: input manifest and expected aggregate values
