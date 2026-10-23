# Header-informed test execution record

Plan: 2026-10-21-header-informed-scientific-tests.md. User authorized implementation with “proceed”.

## Implemented

- Task 1: Full-header TLFB and medication fixtures; weeks/classes/subtypes, repeated visits, partial missingness, positive evidence, paired medication codes/names and dose distractors. Extracted pure functions used by the actual sensitivity script. Fixed positive evidence becoming unknown and inconsistent decimal medication-code representations.
- Task 2: Full-header age/QC/demographic fixtures run preparation scripts 01 and 06. Follow prepared ages through load_modalities(), verify 6-week postnatal age plus GA 38 gives PMA 44, correct run-specific QC, demographic sources, and invariance to competing fields and source-row order.
- Task 3: Full 283-parcel structural/connectivity fixtures; 279 assigned parcels, unprefixed PostDefault, signed network metrics, independent axis permutations and named structural values. Synthetic named CIFTI atlas demonstrates stable effect localization with noncontiguous labels and reordered parcel columns. Regional summary test distinguishes significant thickening, thinning and nonsignificant parcels.
- Task 4: All three scanner schemas; real Python scanner preparation followed by R modality loading/sample selection. Nine identified infants plus unknown identities cannot satisfy the ten-person threshold; ten identified infants do; software versions remain distinct.
- Task 5: Full bundle-statistics schema, competing summaries, reversed subject/bundle patterns, shuffled rows/columns, unrelated metrics. Confirms masked median under the currently coded metric vocabulary.
- Task 6: Full brain-volume schema; unequal tissue components sum to 325, non-tissue and age distractors do not affect it, invalid volumes remain missing, downstream log-volume retains participant identity. Subsequent upstream evidence established Vermis as disjoint tissue: adding 7 must give 332. Production now includes it.
- Task 7: Documentation and user clarification enabled birth-weight units, response ranges and unknown/declined flags, plus item-12 semantics tests. Remaining questionnaire gating and aggregate missingness policies are explicitly deferred.

## Decisions and limits

Ruling: Require atlas label names to match parcel names (after hash removal), instead of retaining an unverified positional fallback. A synthetic test reproduced silent loss/misplacement under positional mapping. Cost: an atlas with generic or differently encoded names now needs a verified name mapping before Figure 4 can render. No actual release atlas mismatch is established by these tests.

Ruling: Preserve existing incomplete-negative and medication-ascertainment policies while correcting conclusive positive substance evidence. Headers do not establish questionnaire skip/unknown coding. Cost: existing classification limitations remain until dictionary review; the new tests do not certify those policies.

Ruling: Use minimal helper extractions for sensitivity covariates and regional network counts, called from the real scripts, instead of running all models to test these transformations. Cost: these tests establish transformation correctness, not full-pipeline statistical calibration.

Ruling: Correct test comparisons to compare participant-aligned values rather than incidental R row names/array attributes. Reordering records changes row names without changing scientific outputs. No expected scientific values were relaxed.

Pending: clarification of the user's control-group construction, medication ascertainment/skip logic, aggregate missingness policy, actual atlas label table, and release diffusion metric units/vocabulary. Age and birth-weight units, birth-weight flags, item-12 meaning and Vermis nonoverlap were resolved using subsequent documentation/user information.

## Verification

- Tests committed before production changes: 3b1299a, dc18345, b97861d.
- Observed RED: missing helper boundaries; parcel mapping returned zero instead of .3 at the specified vertices; positional hemisphere assumption rejected a valid permutation.
- Disposable mutation checks rejected positive-plus-missing aggregation, decimal medication mismatches, adjusted age substituted for candidate age, unknown scanners assigned to a known batch, median substituted for masked median, and unsigned thickening counts. Assertions failed on scientific values or cohort membership, not just exceptions. The first bundle mutation accidentally used pandas' median method; corrected the mutation to column access before accepting its result.
- uv run pytest --cov --cov-report=term-missing -q: 88 passed, 98.89% Python statement coverage.
- uv run make test-r-full: all 32 R test blocks passed after subsequent clarifications, including actual ComBatFamily harmonization.
- Formatting: Black and Air applied to changed files. Unrelated R/05_figures.R formatting restored.
- Existing structural preparation emits pandas fragmentation warnings with full-width fixtures; ComBatFamily emits its existing import-collision warnings. Neither changed asserted results.

## Subsequent documentation and user clarification

The user supplied the public HBCD documentation, resolving several initial deferrals:

- [Age definitions](https://docs.hbcdstudy.org/latest/instruments/agevariables/): tabulated candidate_age is years; adjusted_age is weeks. Updated the competing-field fixture to reflect both units.
- [Pregnancy/infant health](https://docs.hbcdstudy.org/latest/instruments/pregexp/pex/): birth-weight __01 is ounces, __02 pounds. Added a participant-aligned conversion test and extracted the unchanged conversion boundary. Medication dose information is not a reliable ascertainment source.
- [TLFB](https://docs.hbcdstudy.org/latest/instruments/pregexp/su/tlfb/): indices 3–9 are sampled prenatal periods across visits, not consecutive gestational weeks 3–9. ASSIST responses gate administration. Unknown TLFB values therefore cannot establish nonuse without gating information.
- [Missingness documentation](https://docs.hbcdstudy.org/latest/datacuration/phenotypes/): shadow matrices distinguish unknown, declined and logic-skipped responses. This resolves the general blank-field question, but not every instrument's skip rule.
- [BIBSNet volume calculation](https://github.com/DCAN-Labs/BIBSnet/blob/main/src/postbibsnet.py) counts each discrete label separately; the [lookup table](https://github.com/DCAN-Labs/BIBSnet/blob/main/data/FreeSurferColorLUT.txt) assigns Vermis label 172. This establishes nonoverlap. Added a regression requiring 7 units of Vermis to increase total tissue from 325 to 332 and corrected the omitted component. The test failed by exactly 7 before the correction. These upstream sources were inspected on 2026-10-21; actual release provenance remains a separate validation matter.
- User supplied item 12 = “None of these”, coded 0 false/1 true. It must not be added as a positive complication. Clarification of the user's intended control-group coding is pending.
- User supplied birth-weight dropdown ranges: ounces 0–15, pounds 0–20; _dk codes 999 unknown and 777 declined. Tests enforce these ranges and exclude explicit unknown/declined records even when numeric components are present. This conservative contradiction rule prevents unverified measurements entering the analysis. The public page's broader outlier filters do not override the supplied questionnaire coding.

Independent reviewer assessed the implementation through birth-weight extraction as ready within scope, with no critical, important or minor findings. Later Vermis and supplied birth-weight-code changes use their own RED→GREEN regressions. The review did not independently rerun the suites.

Final: Ruling: Retain unresolved questionnaire gating and medication ontology as explicit limits; code normalization does not certify clinical completeness. Cost: incomplete ascertainment may still affect sensitivity samples. The atlas compatibility and statistical-calibration limits above also remain. Item 12 semantics and Vermis nonoverlap are now resolved by subsequent evidence.

Final verification: 88 Python tests passed (98.89% statement coverage); 32 R test blocks passed with ComBatFamily required; git diff --check passed. No tests are skipped for the unresolved definitions; those tests have deliberately not been asserted as implemented. Additional disposable mutations verified that reversed birth-weight units and inversion of the none-of-these response are rejected by value assertions.

## Clarifications incorporated after the initial implementation

The user confirmed that other_complication uses items 1–11 excluding 2 and 3:
1 if any selected item is 1, missing if all selected items are missing, otherwise 0.
Item 12 is excluded. This matches production, so no cohort-definition change was
needed. Expanded the regression to cover all-missing, partial-positive,
partial-negative, all-negative, gestational-hypertension-only and
preeclampsia-only records. The control-group clarification is resolved.

The user clarified that the actual atlas uses Label_n names. The initial
name-only implementation was incompatible with this atlas and has been corrected.
Downloaded the public Myers–Labonte 0.5 atlas and metadata at repository revision
b4403102b96f73e1aec89be8f643cdc6c5172a03. All 283 release-order networks match
parcel_assignments.txt; labels 1–146 are left cortex and 147–283 right cortex,
matching the release prefixes. A fixed numeric-label/name reference now supports
generic labels while retaining invariance to later column permutations. The new
test failed against the name-only implementation before the correction.

Ruling (supersedes name-only requirement): Use the supplied release-order
convention for numeric atlas labels, independently checked at network and
hemisphere level. Cost: those checks cannot independently disambiguate parcels
within one network and hemisphere; the fixed reference preserves the original
release order rather than proving that convention from network memberships alone.

Read Hendrickson et al., Developmental Cognitive Neuroscience 79 (2026) 101706,
DOI 10.1016/j.dcn.2026.101706, from the user-supplied PDF. Sections 2.1.2 and
2.2.4 describe anatomical-label segmentation and resampling; the paper refers to
Supplemental Table 2 for label numbers/names, which is not included in the supplied
main-article PDF. This is consistent with discrete, nonoverlapping regions. The
direct evidence for Vermis specifically remains the upstream label-172 lookup
and per-label voxel counting documented above. The user's Carina check will
confirm the actual release segmentation. No further volume change was made.

Verification after these clarifications: 90 Python tests pass, 98.91% statement
coverage; all 32 R test blocks pass with real ComBatFamily checks.
