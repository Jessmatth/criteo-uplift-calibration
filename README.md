# Calibrating Ad Attribution Against a Randomized Benchmark

Observational advertising measurement, graded against experimental ground truth
on **13.98M rows of real randomized ad-exposure data**.

The dataset is the log of an incrementality test Criteo actually ran: a random
15% of users were withheld from targeting entirely. That randomization makes the
true incremental effect of the campaign computable without assumptions — which
means every observational method can be scored against a number that is not in
dispute.

Written in **R** with the data layer in **SQL** (DuckDB over a 297 MB archive),
reproducible from a clean checkout with one command.

---

## The one-sentence version

Last-touch attribution credits the ad with **23,031** of the 36,711 conversions
among ad-eligible users. The experiment says **10,017** of them were incremental.
Last-touch overstates the ad's contribution by **2.3x**, so the true cost per
incremental conversion is 2.3x what the dashboard reports.

---

## Three numbers that all get called "the effect of the ad"

| Reported as | Effect on conversion | Applies to | Causal? |
|---|---:|---|---|
| Last-touch: exposed vs unexposed | **5.259 pp** | the 3.6% of treated users who saw an ad | no |
| Effect of seeing an ad (CACE, IV) | **3.196 pp** | users whose exposure the campaign changed | yes |
| Effect of running the campaign (ITT) | **0.115 pp** | every user the campaign could reach | yes |

They span a factor of **46**. Which one a team reports decides whether the
campaign looks profitable.

The gap between the first two is pure selection. `exposure` is *post*-randomized:
the ad server chose who got an impression, using signals correlated with buying
intent. Randomized assignment is a valid instrument for exposure, so the Wald
ratio recovers the effect among compliers — and comparing exposed to unexposed
overstates even *that* by 65%.

### The selection, measured rather than assumed

Because assignment was random, the untreated conversion rate of the users the ad
server chose to reach can be backed out of the control arm:

| | Conversion rate |
|---|---:|
| Users the ad server would not reach (never-takers), untreated | 0.119% |
| Users the ad server chose to reach (compliers), untreated | **2.182%** |

Compliers were **18.3x** more likely to convert *before any ad was shown*. That
ratio is the targeting selection, read directly off the experiment. It is the
entire reason last-touch attribution is wrong, and it is normally invisible.

---

## The experiment needed auditing before it could be a benchmark

A randomized benchmark is only a benchmark if the randomization held. Under
clean assignment, standardizing the arm difference to the marginal distribution
of any covariate must leave it unchanged.

It does not. Binning each of the 12 features by value and reweighting to pooled
bin sizes moves the ITT **down** for **9 of 12 features**, by up to 4.6 standard
errors. Standardizing on the full multivariate index moves it from 0.1152 pp to
**0.0843 pp** — a 27% reduction, 8.4 SEs.

Criteo's holdout is not orthogonal to the features (max |SMD| across f0–f11 is
0.049, which at n=14M is roughly 65 SEs from zero). So this study takes the
covariate-standardized quantity as its calibration target and reports the crude
figure alongside it. In conversions, that audit is worth 3,670 conversions of
claimed incremental volume — the difference between believing last-touch
overstates by 1.7x and knowing it overstates by 2.3x.

**This is the study's first finding: the experimental ground truth itself
required observational adjustment before it was usable.**

---

## What gets built

The confounding is manufactured on purpose. An estimator can only be graded on
how much bias it removes if the amount put in is known.

1. **Baseline propensity.** Fit a conversion model on the **control arm only** —
   those users never saw an ad, so their conversion probability is clean
   purchase intent, the signal a real targeting system chases.
2. **Targeting index.** `g(X) = z + 0.60(z² − 1) + 0.45·(f3·f8)`, deliberately
   nonlinear so a linear propensity model is genuinely misspecified.
3. **Strata.** 50 bins on the *value* of `g`.
4. **Selection.** Retain treated users at a rate rising in the stratum, controls
   at a rate falling in it. Strength `alpha ∈ {0, 1, 2, 4}`; intercepts are
   solved so every scenario has the same sample size, making bias comparable
   across the sweep rather than confounded with precision.

Two arms:

- **observed** — targeting driven by features the estimators can see. Bias is
  recoverable in principle.
- **hidden** — targeting also driven by `f0`, which is withheld from every
  estimator. This arm answers *when should an observational readout not be
  trusted*, which method sophistication cannot.

  The loading is calibrated rather than assumed: `f0` is itself 56% predictable
  from the visible features, so hiding it hides nothing by default. At the
  loading used, **30.8%** of the targeting index cannot be reconstructed from
  anything an estimator can see.

Because retention is constant **within** a stratum, sampling leaves the
within-stratum covariate distribution untouched. The benchmark for each
subsample is therefore the stratum-specific randomized differences reweighted to
that subsample's composition — exactly unbiased, with no model in it.

---

## Results

### Every method leaves a large overstatement

Relative bias against the randomized benchmark, conversion outcome, ATE
estimators. **0% means the estimator recovered the experiment.**

**Observed-confounding arm** — every driver of targeting is visible to the estimators:

| Specification | α=0.25 | α=0.5 | α=1 |
|---|---:|---:|---:|
| naive (unadjusted) | 367% | 481% | 520% |
| IPTW [linear] | −80% | 90% | 243% |
| IPTW [flexible] | −311% | 72% | 238% |
| IPTW [oracle] | **7%** | 49% | 240% |
| g-computation [rich] | 30% | 49% | 143% |
| AIPW [linear / rich] | **15%** | 64% | 139% |
| AIPW [flexible / rich] | 30% | 71% | 137% |
| AIPW [oracle / rich] | 37% | 77% | 137% |

**Hidden-confounder arm** — targeting also driven by `f0`, withheld from every estimator:

| Specification | α=0.25 | α=0.5 | α=1 |
|---|---:|---:|---:|
| naive (unadjusted) | 358% | 511% | 588% |
| IPTW [linear] | 34% | 81% | 215% |
| IPTW [oracle] | **4%** | 17% | 265% |
| g-computation [rich] | 55% | 76% | 147% |
| AIPW [linear / rich] | 52% | 81% | 148% |
| AIPW [flexible / rich] | 57% | 85% | 145% |

Three things are worth reading off these tables.

**Adjustment helps enormously and is still not enough.** The naive readout is
360–590% too high. The best realistic specifications pull that to 15–57%. That
is a large reduction and an unusable answer: a 50% overstatement of incremental
conversions is the difference between a campaign that clears its hurdle rate and
one that does not.

**"Share of the gap closed" flatters.** Closing 90% of a 500% gap leaves a 50%
overstatement. This repo leads with relative bias for that reason; the
gap-closed view is in `output/tables/calibration.csv` for anyone who wants it.

**Hiding one targeting driver triples the residual bias.** AIPW [linear/rich]
goes from 15% in the observed arm to 52% in the hidden arm at α=0.25, while the
oracle score barely moves — it knows what the estimators cannot see. That
contrast is the whole point of the arm, and it is what "unmeasured confounding"
costs in a number rather than a caveat.

### Where the ceiling is

Best specification available to a real analyst (oracle excluded). **Selected
after seeing the truth**, so it is an optimistic bound, not a procedure anyone
could follow:

| Arm | α | naive bias | best feasible bias | CI covers truth? |
|---|---:|---:|---:|---|
| observed | 0.25 | 367% | **15%** | yes |
| observed | 0.5 | 481% | 37% | no |
| observed | 1 | 520% | −125% | no |
| hidden | 0.25 | 358% | 34% | no |
| hidden | 0.5 | 511% | 25% | yes |
| hidden | 1 | 588% | 23% | no |

Even picking the winner with knowledge of the answer, only two of the six cells
produce an interval containing the truth. In the other four the reported
uncertainty is saying something false about how much is known.

### Double robustness, measured rather than asserted

The propensity model is misspecified by construction at `ps=linear` — selection
is quadratic plus an interaction, the model is linear. That is exactly the case
AIPW is supposed to survive, and it does. At α=0.25 in the observed arm:

| | relative bias |
|---|---:|
| IPTW with the misspecified score | **−80%** |
| AIPW pairing that same score with a rich outcome model | **15%** |

The protection is one-sided, and the repo says so. Only the propensity model can
be certified correct here, because the conversion surface is real data whose
true form nobody knows. This demonstrates AIPW rescuing a bad propensity model,
not the converse.

IPTW is also the most fragile family by a wide margin. `IPTW [flexible]` at
α=0.25 overshoots to **−311%**: adding the correct functional form made it
worse, because the extra terms produce extreme weights. Being right about the
model and wrong about the answer is a normal outcome, not a pathology.

### Where estimation stops being possible

At α=1 the true propensity score reaches 1.000 and **11.6%** of units fall
outside the truncation bounds — whole strata retain almost no controls. Every
method lands 125–265% off, including the one using the true score.

**No estimator recovers an effect the data cannot identify.** The positivity
guardrail is what says so before the number reaches a slide, and it is the
reason α=1 is in the sweep at all.

### The guardrails, and what actually fails

Share of specifications passing all applicable checks:

| Outcome | base rate | α=0 | α=0.25 | α=0.5 | α=1 |
|---|---:|---:|---:|---:|---:|
| conversion | 0.29% | 0% | 7% | 0% | 0% |
| visit | 4.85% | **63%** | 27% | 7% | 0% |

Failure counts across all conversion specifications: relative bias 118, CI
excludes the benchmark 104, CI too wide 38, balance 28, positivity 24, effective
sample size 0.

**The same machinery is usable for visits and unusable for conversions.** At a
0.29% event rate, 2M rows leave roughly 6,000 conversions carrying the whole
estimate, and at α=0.5 **63%** of conversion specifications fail on interval
width alone — too vague to act on even when the point estimate is close. That is
a precision problem with a different fix (more data, not better adjustment), and
separating it from bias is why both outcomes are analysed.

### The negative control found something real

At α=0 no confounding is induced, so the only thing left to correct is Criteo's
own residual imbalance:

| | observed arm | hidden arm |
|---|---:|---:|
| naive (unadjusted) | 52% | 71% |
| AIPW [oracle / rich] | 30% | 48% |
| IPTW [oracle] | **12%** ✓ covers | **14%** ✓ covers |

IPTW with the true score recovers the benchmark; AIPW does not. This is not a
defect. The oracle score and the benchmark are both defined at **stratum**
resolution, while AIPW's augmentation adjusts at **feature** resolution through
an outcome model whose true form is unknown. Where the underlying randomization
is imperfect, those two resolutions answer slightly different questions.

The pipeline's self-check asserts the IPTW identity, because it is an algebraic
consequence of the construction, and halts the run if it fails.

### Figures

| | |
|---|---|
| `output/figures/relative_bias.png` | relative bias by method and targeting strength, against the 10% tolerance band |
| `output/figures/calibration.png` | every estimate with its interval, against the benchmark |
| `output/figures/balance.png` | worst covariate imbalance before and after weighting |

---

## Two data properties that silently corrupt this analysis

Both were found by assertions in the pipeline, and both are the kind of thing
that produces a plausible-looking wrong answer rather than an error.

**The file is block-ordered by treatment.** The first ~40% of rows are 100%
treated. Row order carries treatment information, so it can never be used as a
tiebreaker.

**Ties are enormous.** 6,684,706 distinct combinations of `f1`–`f11` cover all
13.98M rows, and the single most common combination accounts for **1,699,977
rows (12.2%)**. `ntile` splits such a tie group across strata *in scan order* —
and scan order is treatment-ordered. The first version of this pipeline used
`ntile` and produced a stratum that was **99.96% treated** when randomization
implies 85%. Retention then depended on file position rather than on X, which
destroys both the selection-on-observables premise and the benchmark.

The fix is value-based binning (`cume_dist`), so every tied user lands in the
same stratum and a stratum is a function of X alone. The assertion that caught
it — treated share must be flat across strata — now runs on every execution and
halts the pipeline.

---

## Reproducing

```bash
Rscript run_all.R
```

Downloads a 297 MB archive, builds a DuckDB database, and runs the full sweep.
About 26 minutes on 8 cores after the download, of which `R/06_estimators.R` is
25; every other step is seconds. The suite of 80 unit tests runs at the end.

```bash
N_OBS=500000 Rscript run_all.R     # smaller subsamples, a few minutes
N_FOLDS=2    Rscript run_all.R     # fewer cross-fitting folds
```

Requires R with `DBI`, `duckdb`, `ggplot2`, `testthat`. No compiled dependencies
beyond DuckDB, and no Python.

```
R/00_config.R    paths, sweep grid, guardrail thresholds
R/01_download.R  fetch and verify the archive (MD5, row count)
R/02_build_analytic.R  stage into DuckDB, validate the experiment
R/03_benchmark.R randomized ground truth, CACE, and the randomization audit
R/04_confound.R  build the confounded subsamples and their benchmarks
R/05_attribution.R  last-touch credit vs incremental truth, real data only
R/06_estimators.R   IPTW, matching, AIPW across the sweep
R/07_calibrate.R    grade every estimate against its benchmark
R/08_guardrails.R   SLA pass/fail scorecard
R/09_attribution_synthetic.R  APPENDIX - simulated journeys, quarantined
R/10_figures.R
```

`protocol.md` is the pre-specified study protocol, including the guardrail
thresholds and the amendments log.

---

## Appendix: synthetic multi-touch journeys

Criteo publishes one row per user with a single ad channel and **no touchpoint
sequence**, so first-touch, position-based and time-decay attribution are not
computable from it. `R/09_attribution_synthetic.R` invents a journey for each
converter from a documented generative process in order to run those rules.

The channel mix it reports for organic, email and social is an artifact of that
process and means nothing. The one number that survives is what each rule
credits to **display**, because the display touchpoint is placed from the real
exposure flag and the incremental total comes from the real randomized contrast.

**No number in the tables above comes from that file.**

---

## Data

Criteo Uplift Prediction dataset v2.1, released with:

> Diemert, E., Betlei, A., Renaudin, C., Amini, M-R. (2018). *A Large Scale
> Benchmark for Uplift Modeling.* AdKDD & TargetAd Workshop, KDD 2018.

13,979,592 rows — verified by row count, not taken from the paper. The commonly
cited 25.3M figure refers to the superseded v2.0 release. Licensed
CC-BY-NC-SA 4.0; the archive is downloaded by the pipeline and is not
redistributed here. Code is MIT.
