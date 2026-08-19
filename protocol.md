# Study Protocol

Pre-specified before estimation. Deviations are recorded under Amendments
rather than edited into the text above.

---

## 0. Read this first: what is measured and what is constructed

This study runs on the **Criteo Uplift Prediction dataset v2.1**, the log of a
real randomized incrementality test run by Criteo. Three layers must be kept
apart, and the repository is organised so that they never mix:

| Layer | Status | Where |
|---|---|---|
| Randomized assignment, exposure, visits, conversions, 12 user features | **Real.** Logged by Criteo. | `R/03_benchmark.R`, `R/05_attribution.R` |
| Targeting-induced selection into the observational subsamples | **Constructed by this study.** The selection rule is written down in §5. | `R/04_confound.R`, `R/06_estimators.R` |
| Multi-touch user journeys | **Invented.** Criteo publishes no touchpoint sequence. | `R/09_attribution_synthetic.R` only |

The confounding is deliberately manufactured. That is the point: an estimator
can only be graded on how much bias it removes if the amount of bias put in is
known. Nothing in layer 2 claims to describe how Criteo actually targets.

Layer 3 is an appendix. No number in the README's headline tables comes from it.

---

## 1. Objective

Quantify how far observational advertising measurement lands from experimental
ground truth, and how much of that distance each standard causal-inference
method recovers.

**Primary.** On subsamples where exposure has been made a function of purchase
propensity, estimate the effect of advertising on conversion using last-touch
attribution, IPTW, propensity-score matching and AIPW, and compare each against
the randomized benchmark computed on the same population.

**Secondary.** Repeat on site visits, a 25x more common outcome, to separate
failures of bias correction from failures of precision.

**Operational.** Decide, per method and per severity of targeting, whether the
estimate is fit to inform a budget. §7 fixes the thresholds in advance.

---

## 2. Data source

Criteo Uplift Prediction dataset **v2.1**, released with Diemert, Betlei,
Renaudin and Amini, *A Large Scale Benchmark for Uplift Modeling*, AdKDD 2018.

- **13,979,592 rows**, one per user. Verified by row count on the downloaded
  archive, not taken from the paper. The widely quoted figure of 25.3M refers to
  the superseded v2.0 release.
- 297 MB gzip, MD5 `d2236769ef69e9be52556110102911ec`.
- Licence CC-BY-NC-SA 4.0. The archive is downloaded by `R/01_download.R` and is
  never redistributed in this repository.

Columns: `f0`–`f11`, twelve anonymized dense continuous user features whose
semantics Criteo did not publish; `treatment`, the randomized assignment flag;
`exposure`, whether an ad was actually served; `visit` and `conversion`.

**Known data properties, established in `R/02_build_analytic.R`:**

- 6,684,706 distinct combinations of `f1`–`f11` cover all 13.98M rows. The single
  most common combination accounts for **1,699,977 rows (12.2%)**. Any procedure
  that ranks users must therefore handle very large tie groups explicitly.
- Rows appear in the file in **blocks ordered by treatment** — the first ~40% are
  100% treated. Row order carries treatment information and must never be used
  as a tiebreaker. See §5.

---

## 3. Design

Three estimands are distinguished throughout, because they are routinely
conflated in advertising measurement and they differ here by a factor of 45.

| Estimand | Contrast | Population |
|---|---|---|
| **ITT** | assigned-to-treatment vs held out | every user the campaign could reach |
| **CACE** | Wald ratio, assignment instrumenting exposure | users whose exposure the campaign changed |
| **Naive exposed** | exposed vs unexposed within the treated arm | not a causal contrast |

`exposure` is post-randomization: the ad server chose who received an
impression, using signals correlated with purchase intent. It is therefore never
used as a treatment variable in the calibration study. The analyzed treatment is
always the randomized assignment flag.

**Instrument validity for CACE.** Relevance: assignment moves exposure from 0 to
3.6%. Exclusion: assignment can only act on a user through an ad being shown.
Monotonicity: no user can be exposed while withheld, which holds by construction
because exposure is impossible under control.

---

## 4. Auditing the experiment before using it as a benchmark

A randomized benchmark is only a benchmark if the randomization held. Under
clean assignment, standardizing the arm difference to the marginal distribution
of any covariate must leave it unchanged.

**Pre-specified check** (`R/03_benchmark.R`): bin each of `f0`–`f11` by value,
recompute the arm difference within bins, reweight to pooled bin sizes, and
compare with the crude difference. Systematic movement is evidence of imbalance.

**If imbalance is found**, the covariate-standardized estimate becomes the
calibration target, and the crude difference is reported alongside it as the
cost of not auditing. This is a pre-specified branch, not a post-hoc rescue: the
standardized quantity is the one that is unbiased under the weaker assumption
that assignment is random *conditional on the features*.

---

## 5. Construction of the confounded subsamples

The ad server in this study is a rule this repository writes down.

1. **Baseline propensity.** Fit a conversion model on the **control arm only**
   (1.2M users). Control users never saw an ad, so their conversion probability
   measures pre-existing purchase intent, uncontaminated by any treatment
   effect. This is the signal a real targeting system is built to chase.
2. **Targeting index.** With `z` the standardized linear predictor of that model,
   `g(X) = z + 0.60*(z^2 - 1) + 0.45*(f3_std * f8_std)`.
   The quadratic and interaction terms are there so that a propensity model
   linear in the features is *genuinely misspecified* rather than merely noisy.
3. **Strata.** Cut users into 50 bins on the **value** of `g` using `cume_dist`.
4. **Selection.** Retain treated users at a rate rising in the stratum and
   control users at a rate falling in it. Strength is set by `alpha`; the
   intercepts are solved so that every scenario yields the same sample size,
   which makes bias comparable across the sweep rather than confounded with
   precision.

**Two design choices carry the weight of the whole study:**

**Retention is constant within a stratum.** Every unit in a stratum shares one
retention probability per arm, so sampling leaves the within-stratum covariate
distribution untouched. The benchmark for a subsample is therefore the
stratum-specific randomized differences reweighted to that subsample's stratum
shares — an exactly unbiased quantity with no model in it.

**Strata are cut on the value of `g`, not with `ntile`.** With a tie group of
1.7M rows, `ntile` splits identical users across several strata in scan order,
and scan order in this file is block-ordered by treatment. An earlier version of
this protocol used `ntile` and produced a stratum that was **99.96% treated**;
retention then depended on file position rather than on `X`, which destroys both
the selection-on-observables premise and the benchmark. `cume_dist` places every
tied user in the same stratum, so a stratum is a function of `X` alone. Strata
are consequently unequal in size.

**Assertion.** Because a stratum is a function of `X` and treatment was
randomized, the treated share must be flat across strata. `R/04_confound.R`
stops the pipeline if any stratum deviates by more than 5 percentage points.

### Sweep

The retention ratio between the most and least targeted stratum runs roughly
`exp(2.9 * alpha)`.

| `alpha` | retention ratio | regime |
|---|---:|---|
| 0 | 1.0x | negative control — no selection induced |
| 0.25 | 3.6x | mild targeting, overlap intact |
| 0.5 | 8.7x | moderate targeting, overlap straining |
| 1 | 18.3x | severe targeting, positivity violated |

At `alpha = 0` every estimator, including the unadjusted difference, must land on
the benchmark; failure there indicts the harness, not the estimator.

`alpha = 1` is retained deliberately as a regime where the question cannot be
answered: the true propensity score reaches 0.9996 and 11.6% of units fall
outside the truncation bounds, meaning whole strata retain almost no controls.
Reporting where estimation becomes infeasible is part of the deliverable.

### Arms

- **observed** — the targeting index uses `f1`–`f11` only. Every driver of
  selection is available to the estimators; bias is in principle recoverable.
- **hidden** — the index also loads on `f0`, which is withheld from every
  estimator. This arm exists to answer "when should an observational readout not
  be trusted", which no amount of method sophistication can answer on its own.

  The loading is calibrated, not assumed. `f0` is itself 56% predictable from
  `f1`–`f11`, so withholding it does not automatically hide anything; at a
  loading of 1.0 the two targeting indices correlate at **0.992** and adjusting
  for the visible features very nearly controls the hidden one as well. The
  loading is set to 10, at which **30.8%** of the targeting index cannot be
  reconstructed from the visible features. Returns flatten beyond that because
  of the 56% ceiling.

`f0` is withheld from the estimators in **both** arms, so the arms differ only in
whether the withheld feature drives targeting.

---

## 6. Estimators

All are fit on the visible features `f1`–`f11` only.

| Family | Estimand | Consistent when |
|---|---|---|
| Naive difference in means | ATE | assignment is ignorable unconditionally |
| IPTW, stabilized and truncated to [0.01, 0.99] | ATE | the propensity model is correct |
| 1:1 PS matching on the linear score, with replacement, 0.2 SD caliper | **ATT** | the propensity model is correct |
| g-computation | ATE | the outcome model is correct |
| AIPW, 5-fold cross-fit | ATE | **either** nuisance model is correct |

Matching targets ATT and is graded against the ATT benchmark. Grading it against
an ATE benchmark would charge it for a difference of estimand rather than a
failure of adjustment.

**Propensity specifications.** `linear` (misspecified by construction);
`flexible` (includes the quadratic and interaction basis that generates the
targeting index); `oracle` (the true score, known in closed form because the
selection rule was written by this study). The oracle is not available to a real
analyst and is included only as a control: with the true score, IPTW and AIPW
must recover the benchmark.

**Outcome specifications.** `parsimonious` (linear in the features) and `rich`
(plus the same basis). Neither is claimed to be correctly specified — the
outcome surface is real data and its true form is unknown. Only the propensity
specifications carry a correctness claim.

Crossing three propensity specifications against two outcome specifications is
what turns "AIPW is doubly robust" from a claim into a measurement.

**Inference.** Influence-function standard errors for IPTW and AIPW, treating
the propensity score as known, which is conservative. For matching, a variance
that conditions on the matched sets and prices control reuse; the bootstrap is
**not** offered, because Abadie and Imbens (2008) show it is invalid for
nearest-neighbour matching estimators. Newcombe hybrid-score intervals for
simple risk differences, since the Wald interval is anticonservative at a 0.3%
event rate.

---

## 7. Guardrails

Fixed in `R/00_config.R` before any estimate existed. An estimate that fails any
applicable check is not fit to inform a budget decision, however close its point
estimate looks.

| Check | Threshold | Why |
|---|---|---|
| Relative bias | ≤ 10% of the benchmark | scaled to effect size, so a fixed absolute tolerance cannot flatter a small effect |
| CI covers benchmark | required | an interval that excludes the truth is miscalibrated whatever its centre |
| CI width ratio | ≤ 3x the benchmark's | an estimate can be accurate and still too vague to act on |
| Positivity | ≤ 5% of units outside [0.01, 0.99] | beyond this the question exceeds what the data support |
| Balance | worst absolute SMD ≤ 0.10 | conventional threshold |
| Effective sample size | ≥ 10% of n | balance bought by concentrating weight on a handful of units is not a measurement |

Checks that do not apply to an estimator — a naive difference has no propensity
model — are recorded as `NA` and excluded from the verdict rather than silently
passed.

**A FAIL is usually a finding, not a defect.** The output of interest is not
which estimators pass but where they stop passing.

**One failure stops the pipeline**: if AIPW with the true propensity score at
`alpha = 0` does not cover the benchmark, the construction is wrong and no other
result can be interpreted.

---

## 8. Analysis plan

1. Randomized benchmark and its audit (§4).
2. Attribution credit vs incremental truth, in conversions, on real data only.
3. Confounded subsamples across `alpha` x arm.
4. All estimators on all scenarios, both outcomes.
5. Calibration: bias, relative bias, share of the naive gap closed, coverage.
6. Guardrail scorecard.
7. Appendix: synthetic multi-touch journeys, quarantined.

Primary reporting is the share of the naive gap closed, by estimator, by
targeting strength, split by arm.

---

## 9. Limitations

1. **The features are anonymized.** Nothing here can be interpreted
   substantively. Whether `f0` is a plausible thing for a real ad platform to
   know but an analyst not to have is unanswerable.
2. **The selection mechanism is synthetic.** Real targeting is dynamic, depends
   on auction outcomes and budget pacing, and is not a static function of user
   features. The observed-arm results are therefore an optimistic bound.
3. **Only the propensity model can be correctly specified.** The outcome surface
   is real and its true form is unknown, so the double-robustness demonstration
   is one-sided: it shows AIPW rescuing a bad propensity model with a rich
   outcome model, not the reverse.
4. **The conversion rate is 0.3%.** At `N_OBS = 2,000,000` roughly 6,000
   conversions carry the whole estimate. Precision, not bias, may bind first;
   the visit outcome is analysed for exactly this reason.
5. **Criteo's holdout is not perfectly orthogonal to the features** (§4). The
   benchmark is covariate-standardized in response, which trades an assumption
   of pure randomization for the weaker one of randomization conditional on the
   observed features.
6. **One dataset, one advertiser, one time period.** Nothing here establishes
   how much confounding targeting induces in general.

---

## 10. Amendments

**A1.** Strata were originally cut with `ntile`. Replaced with value-based
`cume_dist` binning after the treated-share assertion in §5 failed at 99.96%.
Cause and consequence documented in §5; the assertion that caught it was added
at the same time and now runs on every execution.

**A2.** The benchmark was originally the crude difference in means. After the
§4 audit found that standardizing on any single feature moves the ITT downward —
by up to 4.6 standard errors, and in the same direction for 9 of 12 features —
the covariate-standardized quantity became the calibration target. The crude
figure is retained in `output/tables/benchmark_randomized.csv` and the gap
between the two is reported as the value of the audit.

**A3.** `visit` was added as a secondary outcome after §9.4 was written, to
separate bias from precision. No thresholds were changed.

**A4.** The sweep originally ran `alpha ∈ {0, 1, 2, 4}`. A first pass showed
that even `alpha = 1` truncated 11.6% of units at the propensity bounds, so
every non-zero cell was a positivity violation and no estimator could recover
the effect in any of them — the guardrails failed everything, which is correct
but uninformative. Rescaled to `{0, 0.25, 0.5, 1}` so the sweep spans feasible
and infeasible regimes rather than only infeasible ones. Verified before
re-running: with the true propensity score, IPTW closes 98% of the naive gap at
`alpha = 0.25`, 90% at `alpha = 0.5` and 54% at `alpha = 1`. No thresholds were
changed.

**A5.** `HIDDEN_LOADING` was raised from 1.0 to 10.0. At 1.0 the observed and
hidden targeting indices correlated at 0.992 and the two arms returned
indistinguishable results, so the hidden arm was not testing what it claimed to.
The loading was chosen by measuring how much of the targeting index remains
unreconstructable from the visible features across candidate values (10.5% at a
loading of 1, 30.8% at 10, 37.3% at 25), not by tuning until the results looked
interesting. No thresholds were changed.
