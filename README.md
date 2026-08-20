# Shared prediction error in Ashokkumar et al. (2026)

Reanalysis code for:

> Wilhelmsen, L. L., Esfandiari, R. & Gollwitzer, A. **Large language models and
> human forecasters make the same prediction errors.** Matters Arising on
> Ashokkumar, A., Hewitt, L., Ghezae, I. & Willer, R. *Large language models can
> predict the results of social science experiments.* Nature (2026).
> https://doi.org/10.1038/s41586-026-10742-x

<!-- Add on acceptance: DOI badge, and the DOI of this repository. -->

## What this is

Ashokkumar et al. show that GPT-4 and a crowd of human forecasters each predict
experimental results about equally well, and propose combining them because they
carry partly independent information. This reanalysis asks about the part not
analysed there: the *errors*. When the model is wrong and the crowd is wrong, are
they wrong on the same experiments, in the same direction?

They are. Across all 1,678 condition comparisons in the primary archive, the two
predictors miss together far more often than independence would produce, and
almost always on the same side. Same-direction errors do not cancel under
averaging, so their agreement is not evidence that a prediction is right.

## Getting the data

The original data are **not redistributed here**. They are available in Code
Ocean capsule **9843791**, released with the original Article. Place these files
in `data/raw/`:

| File | Contents |
|---|---|
| `rct_responses.RDS` | Individual responses of the experiment participants |
| `llm_responses.RDS` | One row per model call, for each of seven models |
| `forecasting_responses.RDS` | One row per forecaster × condition |
| `megastudies.RDS` | The secondary archive of 15 megastudies |

## Running it

```r
# Restore the recorded package versions (recommended)
renv::restore()

# Full pipeline, from the capsule
Rscript run_all.R

# Or: reproduce every reported number from the deposited derived data,
# without obtaining the capsule
Rscript run_all.R --from-derived
```

Roughly 10–20 minutes end to end, dominated by the study-specific null and the
cluster bootstrap. All random components use fixed seeds.

To check the analytical machinery without any data — the corrected-null
integrals against every benchmark reported in the Supplementary Information,
the recalibration helpers against explicit leave-one-study-out loops, and the
cluster-robust sandwich against its expanded form:

```r
Rscript tests/test_functions.R
```

**The pipeline splits where the data stops being redistributable.** Scripts 01
and 02 rebuild the comparison-level datasets from the capsule. Everything
downstream runs from those derived files alone, which *are* deposited here — so
every number in the paper can be reproduced without access to the original data,
and anyone who does obtain the capsule can additionally verify the derivation.

## Layout

```
R/
  00_functions.R          Corrected-null integrals, recalibration helpers,
                          cluster-robust SE, reported-value accumulator
  01_build_comparisons.R  Capsule → observed effects, predictions, recalibration
                          → data/derived/comparisons_primary.csv
  02_build_archive2.R     megastudies.RDS → comparisons_archive2.csv
  03_main_results.R       Miss rates, corrected null, permutation test, study
                          bootstrap, combined predictor, cross-model, experts
  04_robustness.R         Study-specific null, alternative recalibrations and
                          thresholds, study weighting, the authors' sampling
                          scheme, study-level heterogeneity
  05_figures.R            Figures 1 and 2
tests/
  test_functions.R        Checks on the helper functions; needs no data
data/
  raw/                    (empty — see above)
  derived/                Deposited comparison-level datasets
output/
  figures/                F1_same_effects.png, F2_shared_blind_spot.png
  tables/                 Supplementary Tables 1 and 2, reported_values.csv
  sessionInfo.txt
docs/
  reproducible_report.html   Narrated walk-through of the whole analysis
```

## Checking a number in the paper

`output/tables/reported_values.csv` lists every quantity quoted in the
commentary and its Supplementary Information, keyed to where it appears:

| location | label | value |
|---|---|---|
| Main text, para 8 | same-side rate | 0.182 |
| Main text, para 8 | direction share given joint miss | 0.900 |
| SI 1.4 | null direction share | 0.712 |
| … | … | … |

So any figure in the manuscript can be traced to the line of code that produced
it without executing the pipeline.

## Method in brief

Each prediction's error is expressed in units of the observed effect's own
standard error, `z = (prediction − estimate) / s.e.(estimate)`, so that `|z| >
1.96` means the prediction falls outside the experiment's own 95% confidence
interval — a *miss*. Predictions are first recalibrated using the authors' own
leave-one-study-out linear correction, so that the shared error cannot be
dismissed as the overshooting they already document.

The central comparison is against a **corrected null**. Both predictors are
graded against the same estimated effect, and that estimate carries sampling
noise of its own; when it happens to land high, every prediction is "too low"
relative to it. Unrelated predictors would therefore already show some same-side
clustering. The corrected null retains that shared term and makes the
predictions otherwise independent, calibrating each one's spread to reproduce
its observed miss rate. Whatever the data show beyond that benchmark is genuine
dependence between the predictions.

Because comparisons within a study share participants and often a reference
condition, all uncertainty intervals come from a cluster bootstrap that
resamples studies rather than comparisons.

## Reproducibility notes

- Paths are resolved with `here::here()` from the project root. There is no
  `setwd()` anywhere; open the `.Rproj` or run from the root.
- Condition names are recoded to ASCII before regression and string comparisons
  use a fixed ICU collator, so results are identical under any locale. Without
  this, R silently mangles typographic quotes in coefficient names on non-UTF-8
  systems.
- `renv.lock` records the exact package versions used.
- `output/sessionInfo.txt` is written on every run.

## Citation

<!-- Replace with the Zenodo DOI on deposit. -->

```
Wilhelmsen, L. L., Esfandiari, R. & Gollwitzer, A. (2026).
Code for "Large language models and human forecasters make the same prediction
errors". Zenodo. https://doi.org/XXXXX
```

## Licence

Code is released under the MIT Licence (see `LICENSE`). The derived
comparison-level datasets in `data/derived/` are released under CC BY 4.0.

The original data remain under the licence of Code Ocean capsule 9843791 and are
not redistributed here.

## Contact

Live Leonhardsen Wilhelmsen — live.l.wilhelmsen@bi.no
Anton Gollwitzer — anton.gollwitzer@bi.no

Center for Democracy and Information Integrity, BI Norwegian Business School
