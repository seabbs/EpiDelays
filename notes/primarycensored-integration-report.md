# primarycensored integration report

Branch: `primarycensored-integration`.
Upstream issue: [oswaldogressani/EpiDelays#1](https://github.com/oswaldogressani/EpiDelays/issues/1).

## 1. Summary

`kerlikelihood()`'s doubly interval-censored branch is, up to a constant that drops out of optimisation, identical to `primarycensored::pprimarycensored` with a uniform primary event window.
Replacing the per-row `stats::integrate()` calls with `pcens_cdf()` gives exact results for gamma, log-normal and Weibull via analytical solutions already shipped upstream, and keeps numerical integration as the fallback for gaussian and skew-normal.
The most substantive port candidate from EpiDelays into primarycensored is the non-parametric delay likelihood (see epinowcast/primarycensored#218); the pure-R skew-normal helpers are a weaker candidate and probably not worth landing upstream.
Benchmarks should confirm that the analytical families see the largest speedup.

## 2. Math equivalence

Take one row of the doubly interval-censored data with primary bounds `[x1l, x1r]` and secondary bounds `[x2l, x2r]`.
Write `a = x1r - x1l` for the primary window width.
The inner integral evaluated at `R/kerlikelihood.R:120` (gaussian), `:174` (skew-normal), `:225` (gamma), `:276` (log-normal), `:327` (weibull) is:

```
I = integrate(h, x1l, x1r),  h(x1) = F(x2r - x1) - F(x2l - x1)
```

Substitute `p = x1 - x1l`, so `p in [0, a]` and `x1 = x1l + p`:

```
I = integrate_0^a [ F((x2r - x1l) - p) - F((x2l - x1l) - p) ] dp
```

With uniform primary density `f_primary(p) = 1/a` on `[0, a]`, primarycensored defines:

```
F_cens(q; a) = integrate_0^a F(q - p) * (1/a) dp = (1/a) * integrate_0^a F(q - p) dp
```

So `(1/a) * I = F_cens(x2r - x1l; a) - F_cens(x2l - x1l; a)`, i.e.

```
log(I) - log(a) = log( F_cens(x2r - x1l; a) - F_cens(x2l - x1l; a) )
```

The left-hand side is exactly what `kerlikelihood` accumulates: `logint - logdx1[i]` with `logdx1 <- log(x$x1r - x$x1l)` prepared at `R/kerlikelihood.R:78`.
The right-hand side is exactly `log( pprimarycensored(x2r - x1l, pwindow = a) - pprimarycensored(x2l - x1l, pwindow = a) )` with `L = 0`, `D = Inf`, `dprimary = dunif`.
This matches the numerical check already carried out for the gaussian family.

## 3. Family coverage matrix

S3 methods present in `primarycensored/R/pcens_cdf.R` (grepped directly):
`pcens_cdf.default`, `pcens_cdf.pcens_pgamma_dunif`, `pcens_cdf.pcens_plnorm_dunif`, `pcens_cdf.pcens_pweibull_dunif`.

| EpiDelays family | EpiDelays implementation | primarycensored analytical? | After port |
|---|---|---|---|
| `gamma` | `stats::pgamma` + `stats::integrate` per row | Yes: `pcens_cdf.pcens_pgamma_dunif` | Analytical, direct swap |
| `lognormal` | `stats::plnorm` + `stats::integrate` per row | Yes: `pcens_cdf.pcens_plnorm_dunif` | Analytical, direct swap |
| `weibull` | `stats::pweibull` + `stats::integrate` per row | Yes: `pcens_cdf.pcens_pweibull_dunif` (falls back to numeric when `pwindow > 3`) | Analytical for `pwindow <= 3`, numeric otherwise |
| `gaussian` | `stats::pnorm` + `stats::integrate` per row | No | Numeric via `pcens_cdf.default`, still benefits from the shared machinery |
| `skewnorm` | `pskewnorm` (pure-R Owen's T) + `stats::integrate` per row | No; primarycensored has no skew-normal class at all | Needs `add_name_attribute` + numeric default; see port candidates below |

Note: primarycensored's Weibull method delegates to numeric integration when `pwindow > 3` (`R/pcens_cdf.R:284`), so any EpiDelays data with a three-day-or-wider primary window keeps numerical integration for Weibull.
This is a known caveat rather than a regression.

## 4. Port candidates

### 4.1 Skew-normal support in primarycensored

EpiDelays ships its own pure-R skew-normal CDF/PDF/quantile functions at `R/pskewnorm.R`, `R/dskewnorm.R`, `R/qskewnorm.R`, using Owen's T function to avoid a dependency on `sn`.
primarycensored has no skew-normal coverage.
The skew-normal distribution does not have a closed-form primary-censored CDF under a uniform primary window in general, so any skew-normal S3 method on `pcens_cdf` would still go through numerical integration.
That means the only thing a port would actually add is a named `pskewnorm` that primarycensored users can pass into `pprimarycensored(..., pdist = pskewnorm)` without pulling in `sn`.
This is genuinely minor.
The honest recommendation: do not port the helpers into primarycensored core; instead, note in primarycensored docs that `sn::psn` works out of the box through the numeric default, and that users who want an `sn`-free path can lift EpiDelays' Owen's T implementation.
If upstream wants a small vignette showing skew-normal through primarycensored, that is a lighter lift than adding a new dependency surface.

### 4.2 Non-parametric likelihood

EpiDelays has `nonparfit()` (`R/nonparfit.R`), a non-parametric estimator based on uniform mixtures (Gressani and Hens, 2025) that operates on single interval-censored data without a parametric family.
primarycensored issue [epinowcast/primarycensored#218](https://github.com/epinowcast/primarycensored/issues/218) already asks for non-parametric delay likelihoods; seabbs raised the cross-link in the EpiDelays issue thread.
The sketch: `pcens_cdf`'s dispatch is parameter-free at the interface — it takes a `primarycensored` object and returns CDF values.
A non-parametric delay CDF can be represented as a step or piecewise-linear function and passed in as `pdist`, with the numerical default handling the primary-censoring integral.
What primarycensored currently lacks is (a) a principled way to tag such a `pdist` so `pcens_cdf` can dispatch to a dedicated method when the delay is piecewise-constant on intervals that coincide with the primary window, and (b) a fitting entry point that updates the non-parametric CDF under the censoring model.
The uniform-mixture construction in `nonparfit` is a natural match because uniform mixtures compose cleanly with a uniform primary density: the convolution stays piecewise polynomial.
A port would add a `pdist` builder that exposes a non-parametric CDF tagged for primarycensored and a fitter that iterates between updating mixture weights and applying the primary-censored likelihood.
This is the most substantive upstream contribution available from EpiDelays.

### 4.3 Truncation in upstream EpiDelays

primarycensored already handles `L` and `D` via `.normalise_cdf` in `primarycensored/R/pprimarycensored.R`.
EpiDelays does not currently expose truncation bounds at all.
The proposed API for adding truncation to EpiDelays lives in a separate design note at `notes/truncation-design.md`.
Once EpiDelays rides on top of `pprimarycensored`, truncation comes for free by forwarding `L`/`D` through `kerlikelihood`'s public entry points (`parfitml`, `parfitmom`).

## 5. Numerical stability and performance notes

Measured on the issue #1 reprex (integer-day data, `n = 100`, one primary window of width 1, one secondary window of width 1), routing through `primarycensored::dprimarycensored` on the full `x2l - x1l` vector gives:

- gamma: ~5.6x faster than the current `kerlikelihood` per-row `stats::integrate` path.
- gaussian: ~2.7x faster.

See the benchmark log at [oswaldogressani/EpiDelays#1 (comment 4235984297)](https://github.com/oswaldogressani/EpiDelays/issues/1#issuecomment-4235984297).

The two speedups come from different mechanisms:

- Gamma wins because `pcens_cdf.pcens_pgamma_dunif` (`primarycensored/R/pcens_cdf.R`) ships a closed-form solution for the primary-censored CDF under a uniform primary window.
  No numerical integration happens at all; the per-row `stats::integrate` in `kerlikelihood` is replaced with arithmetic on `pgamma` evaluations.
  Log-normal should behave similarly once benchmarked, and Weibull with `pwindow <= 3` will too.
- Gaussian wins because `dprimarycensored` dedupes quantile evaluations.
  Internally it computes `unique_points <- sort(unique(c(x, x + swindow)))` and calls `pprimarycensored` once on that set, caches the CDF values in a lookup table, and differences them per row (`primarycensored/R/dprimarycensored.R:110-145`).
  For integer-day data the set of unique points is tiny relative to `n`, so even though the default gaussian path still uses `stats::integrate` under the hood, it runs on far fewer quantiles than the per-row loop does.

**Lesson from the first swap iteration.** An earlier version of the `primarycensored` engine in `R/kerlikelihood.R` called `pprimarycensored(c(ql, qr), ...)` once per row.
This is correct but it defeats the lookup-table dedup inside `dprimarycensored`, and leaves roughly 4–6x on the table on integer-day data.
The idiomatic path is a single vectorised `dprimarycensored(x = x2l - x1l, pdist, pwindow, swindow, log = TRUE)` call (grouped by unique `(pwindow, swindow)` when those vary across rows).
The engine is being rewritten in parallel to use this pattern.

**Weibull papercut still stands.**
`pcens_cdf.pcens_pweibull_dunif` silently falls through to the numeric default when `pwindow > 3` (`primarycensored/R/pcens_cdf.R:284`).
For the common case `pwindow == 1` it takes the analytical path and behaves like gamma and log-normal.
It is worth raising the silent fallback upstream, either as a warning or a documented threshold, since users running Weibull with wider primary windows will see a performance cliff they cannot predict from the API surface.

Stability should be at least as good as the current implementation for the analytical families since they avoid cancellation errors from `integrate` near tails; the default numeric path matches the existing behaviour modulo integration tolerances.

A follow-up issue should benchmark `parfitml` on all five families across representative sample sizes once the rewrite lands, so the speedup table can be reported end-to-end rather than at the per-call level.

### End-to-end parfitml benchmark

End-to-end wall-clock from a single `parfitml()` fit per cell on doubly-interval-censored data, with `pwindow = swindow = 1` and `Bboot = 100`.
Data are simulated via `primarycensored::rprimarycensored()` using the same ground-truth parameters as `tests/testthat/test-kerlikelihood-primarycensored.R`; for the non-negative-support families (`gamma`, `lognormal`, `weibull`) draws with `delay < pwindow` are rejected so `parfitml`'s non-negative sample-space check passes on all rows.
Reproduce with `Rscript scripts/benchmark-parfitml.R` from the repo root.

| family    |    n | Bboot | elapsed (s) | converged |
|-----------|-----:|------:|------------:|-----------|
| gaussian  |   50 |   100 |        1.91 | TRUE      |
| gaussian  |  200 |   100 |        3.66 | TRUE      |
| gaussian  | 1000 |   100 |       10.48 | TRUE      |
| gamma     |   50 |   100 |        2.18 | TRUE      |
| gamma     |  200 |   100 |        3.55 | TRUE      |
| gamma     | 1000 |   100 |       10.80 | TRUE      |
| lognormal |   50 |   100 |        1.93 | TRUE      |
| lognormal |  200 |   100 |        3.28 | TRUE      |
| lognormal | 1000 |   100 |       10.50 | TRUE      |
| weibull   |   50 |   100 |        1.90 | TRUE      |
| weibull   |  200 |   100 |        3.41 | TRUE      |
| weibull   | 1000 |   100 |       10.81 | TRUE      |
| skewnorm  |   50 |   100 |       38.17 | TRUE      |
| skewnorm  |  200 |   100 |       29.08 | TRUE      |
| skewnorm  | 1000 |   100 |       41.82 | TRUE      |

The four numeric families land within a narrow band of each other at each sample size, which confirms that the per-row `dprimarycensored` dispatch inside `build_pc_loglik` is no longer the bottleneck; the analytical `pcens_cdf` fast path for gamma, log-normal and Weibull pulls them into the same regime as gaussian despite gaussian still routing through numerical integration.
Scaling is sub-linear in `n` relative to a pure `n × Bboot` model: going from `n = 50` to `n = 1000` (a 20x jump in data) multiplies wall-clock by about 5 to 6, so the cost is dominated by a mix of `optim` iterations and bootstrap resamples rather than by the per-row likelihood evaluation, consistent with `dprimarycensored`'s unique-quantile dedup flattening the per-call cost on integer-day data.
Skew-normal now converges on all three cells after two fixes: `kerlikelihood()` wraps `pskewnorm` with a clamp to `[0, 1]` plus a `cummax` on sorted inputs so `primarycensored::check_pdist` does not reject the CDF when Owen's T saturates mid-optim, and `parfitmom()` uses a signed cube root (plus a finiteness fallback) so small-sample negative skew does not seed the optimiser with NaN.
Skew-normal is slower than the other families because `pskewnorm` does per-quantile Owen's T evaluations and cannot share the analytical `pcens_cdf` fast path; per-call cost dominates the benchmark cell.

Benchmark machine: R 4.5.0 (2025-04-11), `aarch64-apple-darwin20`.
Total wall-clock for all 15 cells: 173.8 s.

## 6. Open questions for upstream

- Would upstream accept a PR that adds an `engine = c("integrate", "primarycensored")` switch to `kerlikelihood`, `parfitml` and `parfitmom`, defaulting to the existing `integrate` code for one release and flipping the default once benchmarks and a CRAN round-trip have landed?
- Is there appetite for exposing `L` and `D` truncation bounds on `parfitml`/`parfitmom` as part of the same PR, or should that ride in a follow-up once the engine switch is merged?
- For non-parametric support, would upstream prefer a thin `nonparfit` re-skin that delegates to primarycensored once epinowcast/primarycensored#218 lands, or to keep `nonparfit` self-contained and contribute the uniform-mixture construction as a reference implementation into primarycensored instead?
- Should `dprimarycensored` accept vector `pwindow` and `swindow` so that clients with row-varying windows don't have to manually group by the unique `(pwindow, swindow)` pairs and re-invoke it per group?
  At present any client with mixed windows either loops per row (losing the dedup speedup) or reimplements the group-by dance externally; pushing the grouping into `dprimarycensored` would restore the analytical/dedup fast path without forcing that boilerplate on callers.
  Worth raising as a concrete primarycensored feature request.
