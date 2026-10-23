# Header-informed Scientific Validity Test Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Use the supplied release headers to expose silent errors in phenotype definitions, participant selection, measurement selection, network calculations, and anatomical attribution.

**Architecture:** Build synthetic tables with the actual columns in `hbcd_headers.json`, filling competing fields with deliberately different values. Exercise production readers and transformations, then assert independently calculated values and participant identities. Retain the existing small mathematical fixtures; these new tests cover ambiguities those fixtures cannot express.

**Tech Stack:** Python pytest, native R testthat, existing NumPy/pandas/nibabel and R dependencies; uv for local commands.

**Spec:** User request of 2026-10-21: plan additional tests using the supplied headers, focusing on scientific validity rather than operational behavior. This document records that scope and the evidence supporting the proposed tests.

Implementation status and subsequent documentation/dictionary findings are recorded in [the execution record](2026-10-21-header-tests-results.md).

## Global constraints

- Implementation authorized by the user's subsequent “proceed”; test-first fixes remain required.
- Use invented participant records, never reconstruct or request actual observations.
- Test silent changes in scientific outputs, not missing-file errors, CLI behavior, exact header counts, or formatting.
- Use pytest for Python and testthat for R, functions rather than test classes, and fixtures for reusable objects.
- Run local commands through `uv run`; do not add imports or code to `__init__.py`.
- Before any production fix, demonstrate a scientifically justified failing test and commit the test. Do not adjust expected results simply to match current behavior.
- Existing correct behavior can be tested with a deliberate mutation in a disposable copy, followed by the unchanged production implementation passing.
- Header names establish structure, not units, valid response codes, skip logic, or the intended estimand. Mark unresolved scientific definitions explicitly rather than encoding guesses as truth.
- Prefer small pure-function extractions if top-level R analysis scripts prevent testing their transformations. Avoid evaluating fragments selected by source-line numbers or substituting fake harmonization.

## Evidence from the supplied headers

All 18 reported file types are present. Scanner metadata has three schema variants; other file types have one each in this report.

- Thickness and sulcal depth each contain 283 parcel columns. Connectivity has `Node` followed by the same 283 parcels, and its row labels match those columns. After removing hash suffixes, all three parcel orders agree. This does **not** verify the separate CIFTI atlas label mapping.
- Parcel 1 is `PostDefault`, with no hemisphere/index prefix. Four parcels are `L_None_1`, `L_None_2`, `R_None_3`, and `R_None_4` after hash removal.
- Diffusion and structural QC both contain `candidate_age` and `adjusted_age`. Diffusion QC contains raw, masked-raw, T1, and masked-T1 neighboring-correlation measures; structural QC has summary scores and multiple view-specific scores/reviewer counts.
- Demographics includes maternal age at Visit 1 and delivery, several maternal/child race definitions, and recruitment site. The separate visit table also supplies site.
- TLFB includes aggregate substance-week fields, subtype fields, and more substance classes than the four currently modeled.
- Pregnancy has `compl_001___1` through `___12`; preparation currently selects only 1–11. Medication and infant-health tables include companion `_dk` fields.
- BIBSNet includes `Vermis`, CSF, ventricles, choroid plexus, and `Unknown`, in addition to the tissue fields currently summed. `Vermis` is not in the current sum.
- One scanner schema has `filename`, acquisition/QC fields, and site but none of the four `scanner_*` identity fields. The other two contain those fields.
- Bundle statistics contains `mean`, `median`, `masked_mean`, `masked_median`, weighted alternatives, and acquisition/provenance fields.

## Review focus

1. A positive substance report combined with missing reports must not silently disappear from a logically defined “any use” group (Task 1).
2. Alternative age/QC columns must not change postmenstrual age or covariate adjustment (Task 2).
3. Same-hemisphere parcel permutations must not silently relocate an effect (Task 3).
4. Missing scanner identity must not create a fictitious common batch or inflate the batch-size threshold (Task 4).
5. Competing summary measures, missing responses, and unconfirmed questionnaire definitions must not be treated as interchangeable (Tasks 5–7).

## Shared fixture approach

Add `tests/header_fixtures.py` with `release_columns(kind, variant=0)` and `release_frame(kind, rows, variant=0)`. Read the committed JSON directly; preserve column names/order, create all columns, and populate only explicitly supplied synthetic values. Expected outputs must not be computed using production selectors.

For R, add `tests/testthat/helper-header-fixtures.R` with `release_columns(kind, variant=1L)` and `release_frame(kind, rows, variant=1L)`. Add `tests/export_header_columns.py`, a standard-library utility invoked with `uv run python`, to export `(kind, variant, column_index, column_name)` records from the committed JSON to a temporary TSV once per R test session. The R helper reads that TSV without name repair and orders by `column_index`. This avoids a new R JSON dependency and manual duplication of hundreds of headers. Keep all fixture data under existing temporary-release isolation.

## Task 1 — Substance-use and medication phenotypes (highest priority)

**Files:** new `tests/testthat/test-sensitivity-phenotypes.R`, shared R header helper; production boundary is covariate construction in `R/03_sensitivity.R`.

**Interfaces:** synthetic full-header TLFB/medication tables in; participant-indexed `tlfb_*`, `any_substance`, and `antihypertensive` values out. If necessary, extract these transformations into `R/sensitivity_covariates.R` and call them from the real script, keeping model fitting outside the helper.

- [ ] Test `prenatal_use_ignores_preconception_and_unselected_fields`: an aggregate nicotine value only in week 02 must not classify prenatal nicotine exposure; the same value in week 03 or 09 must. Populate actual subtype and other-class columns with distracting values. These tests pin the currently stated four-class/week-3–9 definition, not its clinical correctness.
- [ ] Test `any_use_retains_positive_evidence_with_partial_missingness`: four-class values `(1, NA, 0, NA)` must give 1; `(0, 0, 0, 0)` gives 0; all missing gives NA. Current `rowSums` without missing-value handling is expected to fail the first case. A `(0, NA, 0, 0)` rule needs an explicit missingness policy before adding its assertion.
- [ ] Test grouping across repeated visits: positive evidence on either visit remains positive; reordering visit rows leaves the participant-level exposure unchanged; an all-missing participant is not classified as a known nonuser.
- [ ] Test medication code/name pairing using an invented code `900001` associated with `labetalol` in the actual `...003__01`/`...003__01_name` pair. Its occurrence in a Visit 2 `...047__01` slot identifies treatment. The same numeral in a dose/timing field does not. Code strings `900001` and `900001.0` should identify the same drug. Use literal known codes only to test code handling, not to independently validate drug ontology.
- [ ] Assert the resulting within-HDP comparison membership, not just Boolean helper outputs: unknown medication ascertainment must follow the agreed policy rather than silently entering “untreated.” Gate `_dk` and questionnaire-skip expectations on the dictionary.
- [ ] Run `uv run make test-r-full`; demonstrate targeted failures or mutations, commit tests, then implement only necessary agreed fixes and rerun.

**Scientific consequence:** changes substance-adjusted estimates, exclusion analyses, and the treated-versus-untreated HDP contrast. The headers also show that the label “any prenatal substance use” currently covers only four classes; confirm whether that restricted definition is intended before expanding it.

## Task 2 — Competing ages, QC measures, and demographic definitions

**Files:** extend `tests/testthat/helper-preparation.R`, `test-preparation.R`, `helper-covariates.R`, and `test-covariates.R` with full-header fixtures.

**Interfaces:** actual preparation scripts 01 and 06, followed by `load_modalities()`, consume temporary release tables and emit model inputs.

- [ ] Give `candidate_age` the synthetic value `42 / 365.25` and `adjusted_age` the value `14 / 365.25`; with delivery GA 38, assert postnatal age 6 weeks and PMA 44 under the code's current years-to-weeks contract. Changing adjusted age alone must not change either result. Confirm actual age units separately before claiming release-unit validation.
- [ ] Set raw neighboring correlation to .81, masked raw to .91, T1 to .71, and masked T1 to .61. Assert model inputs contain .81. Set BrainSwipes summary QC to 2 while view scores and review counts differ; assert the summary enters the model.
- [ ] Supply two structural runs with different ages/QC values; verify each model row receives values from its selected run. Compare participant IDs and values after source-row reordering.
- [ ] Set maternal delivery age 30 versus Visit 1 age 29; give maternal AOU race, maternal screening race, child race, and the two site sources distinguishable values. Assert the exact source fields currently specified by preparation, while retaining existing earliest-nonmissing demographic/latest-visit-site rules.
- [ ] Vary unused competing columns and require identical model-input tables. Run and commit tests before any production fix.

**Scientific consequence:** prevents double age adjustment, wrong quality adjustment, and accidental changes in the demographic adjustment set.

## Task 3 — Actual parcel names, networks, and effect localization

**Files:** new `tests/test_header_parcels.py`; extend `tests/test_figures.py`; production boundaries are `prepare/04_structural.py`, `prepare/05_functional.py`, `python/figure4_map.py`, and regional network summaries.

- [ ] Build synthetic tables using the full 283-name schema, retaining the real `PostDefault` name, hashes, and four `None` parcels. Use a hand-calculated active subgraph on named parcels and zero-strength remaining nodes. Assert network assignments by name, exclude unassigned parcels from functional metrics, and retain the unprefixed `PostDefault` parcel in its intended network. With complete synthetic coverage, expect 279 retained assigned parcels; average node participation over those 279 nodes, not only the active subgraph. Do not hard-code the real analysis's 276 retained parcels into this complete-coverage fixture.
- [ ] Shuffle connectivity rows and columns independently and permute structural columns; require unchanged named values, signed network statistics, and effect locations. The reference parcel-to-label mapping must be held fixed independently of the permuted input.
- [ ] Give two actual left-hemisphere parcel names effects +.3 and -.2. With a synthetic CIFTI label table keyed to those names, assert the correct vertex receives +.3 and the negative-effect region stays absent from the thicker-cortex map. A same-hemisphere swap must not pass unnoticed.
- [ ] Test thicker/thinner network counts using these names: one significant positive, one significant negative, and one nonsignificant parcel should contribute counts 1, 1, and total 3 within the named network.
- [ ] Run `uv run pytest tests/test_header_parcels.py tests/test_figures.py -q`; prove sensitivity to positional mapping and broken network-name parsing before implementing fixes.

**Limit:** the supplied TSV schemas agree with each other. Do not claim an actual atlas mismatch. Verification of the real figure's anatomical assignment requires the CIFTI numeric-label-to-name mapping, which is not in this header report. Synthetic tests can establish invariance but cannot supply that missing reference. Regional handling of unassigned parcels also needs an explicit scientific policy before asserting a change in the regional testing family.

## Task 4 — Scanner schema variants and analytic membership

**Files:** new `tests/test_header_scanners.py`; extend `tests/testthat/test-common.R` for the downstream batch consequence.

- [ ] Generate files with each of the three real scanner schemas. For the two identity-bearing variants, identical synthetic serial/software values must yield the same batch despite differing unrelated QC columns.
- [ ] For the variant lacking identity fields, never substitute site, age, or QC as scanner identity. Assert those participants remain unknown and do not join a fabricated scanner batch.
- [ ] Follow the emitted metadata through modality selection: nine identified participants plus one unidentified participant must not meet the ten-person batch threshold; ten identified participants plus one unidentified participant retain exactly the ten identified IDs.
- [ ] Test distinct software versions on the same serial number remain distinct where the current batch definition requires this. Do not infer that two acquisitions are the same scanner merely because their site agrees.
- [ ] Verify cohort membership and batch assignments, not just parsing or absence of exceptions; run Python tests plus `uv run make test-r-full`.

**Scientific consequence:** harmonization groups and exclusion of scanner-confounded observations. Headers establish that identity can be unavailable, not how often it is unavailable in the actual analytic cohort.

## Task 5 — Bundle-summary identity

**Files:** extend `tests/test_prepare.py` and the shared Python header helper; production boundary `prepare/03_diffusion.py`.

- [ ] Populate all 23 real bundle-statistics columns. Give `masked_median=.2`, `median=.7`, `mean=.8`, and weighted summaries=.9 for one FA row. Assert the emitted named FA bundle value is .2.
- [ ] Give MD and MKT distinct values and shuffle rows/columns; require correct participant, bundle, and metric association. Use unrelated metrics as distractors so selecting the first numeric field cannot pass.
- [ ] Include two synthetic subjects with opposite bundle patterns and verify no cross-participant pairing after pivoting.
- [ ] Run `uv run pytest tests/test_prepare.py -q`; mutation-check selection of `median` instead of `masked_median`.

**Limit:** headers do not reveal the actual `variable_name` vocabulary or the units of FA/MD/MKT. Confirm those in documentation before adding tests that purport to validate tensor-versus-DKI model provenance. Do not fabricate that vocabulary from the production filter as an independent oracle.

## Task 6 — Total brain volume and competing anatomical components

**Files:** extend `tests/testthat/helper-covariates.R` and `test-covariates.R` with all 40 volume columns.

- [ ] Use unequal left/right and structure-specific synthetic volumes. Assert a hand-summed tissue total, not the current uniform-value fixture alone.
- [ ] Vary CSF, ventricles, choroid plexus, `Unknown`, and age fields by large amounts without changing the tissue total. Preserve existing missing/negative/empty-hemisphere tests.
- [ ] Confirm whether `Vermis` is a disjoint cerebellar tissue component that should be included in the intended total. If yes, adding 7 units of Vermis must add exactly 7 units to total volume; if it overlaps an included component, require no double counting. Do not lock in the current omission without that decision.
- [ ] Assert the downstream `log_tbv` value against the hand-summed total and check its participant association.
- [ ] Run `uv run make test-r-full`; commit the agreed volume-definition regression before changing the sum.

## Task 7 — Questionnaire definitions requiring dictionary confirmation

**Files:** extend `tests/testthat/test-covariates.R` and `test-sensitivity-phenotypes.R` after definitions are confirmed.

- [ ] Obtain the label for `pex_bm_healthv2_preg__compl_001___12`. If it is a relevant complication, a positive item 12 alone must affect the composite; if it is a “none/other response” control, test its documented treatment and contradictions. Its existence alone is not evidence it should be added.
- [ ] Confirm units/meanings of infant weight `...001__01` and `...001__02`, and the precedence of `...001__dk`. Under the currently coded ounces/pounds mapping, 7 lb 8 oz must yield 3401.94 g. A reversed assignment must fail. Do not claim this is the correct release mapping until verified.
- [ ] Confirm `_dk` codes, nulliparity skip rules, and labor-onset/cesarean category meanings. Test known positive, known negative, unknown, and contradictory value/flag combinations, asserting downstream inclusion and covariate values.
- [ ] Establish one missingness policy for aggregate exposure variables: positive evidence with partial missingness, negative evidence with partial missingness, and entirely unknown records. Then test that policy consistently across diabetes, complications, substance use, aspirin, and medication treatment where scientifically applicable.

## Execution order and completion criteria

Implement Tasks 1–5 first where expectations are already specified, while collecting the definition decisions needed for Tasks 6–7. Keep tests pending on unresolved meanings out of the passing-suite claim; do not silently skip them and call the plan complete.

For each task: write tests, observe justified failure or a targeted mutation failure, commit tests, make only the required production change, and verify green. Existing tests that already cover the behavior should be extended rather than duplicated. A failing scientific test is a finding to report, not a reason to weaken the test.

Final checks: `uv run pytest --cov --cov-report=term-missing -q`, `uv run make test-r-full`, and `uv run git diff --check`. Add targeted disposable mutations for competing-column selection, partial-missingness aggregation, same-hemisphere parcel swaps, and scanner-batch contamination. No coverage-percentage increase alone constitutes evidence of scientific validity.

This plan deliberately excludes shell behavior, missing-file diagnostics, JSON formatting, and general dependency/CLI tests. It also does not substitute for full-pipeline calibration of harmonization or bootstrap inference, which the headers cannot establish.
