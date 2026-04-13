# Design: Truncation support via primarycensored integration

Issue: oswaldogressani/EpiDelays#1.
Scope: add left/right truncation (`L`, `D`) to `parfitml()` and the underlying kerlikelihood factories by routing the double-censored likelihood through `primarycensored::pprimarycensored()`.
This note is design-only; no code is changed.

## 1. API proposal

### Public signature

```r
parfitml(
  x,
  family,
  Bboot  = 1000,
  pgbar  = TRUE,
  L      = 0,
  D      = Inf
)
```

`L` and `D` are scalars applied to every row.
Row-varying truncation is a later extension (see open questions).
Validation follows `primarycensored::.check_truncation_bounds`: `L >= 0`, `L < D`.
Additional EpiDelays-side checks (run before fitting):

- `nc == 4` branch: all rows must satisfy `x2l - x1r >= L` and `x2r - x1l <= D`.
  A row whose observed delay interval cannot fall inside `[L, D]` contributes zero mass and yields `-Inf` log-likelihood; reject at input time with a clear message.
- `nc == 2` branch: `xl >= L` and `xr <= D`.
- Non-negative-support families (gamma, weibull, lognormal) already require `xl >= 0`; `L >= 0` is consistent.

`parfitmom()` gains the same `L`, `D` arguments so that `parfitml()` can forward them (see 4).
`kerlikelihood()` gains `L`, `D` arguments that are captured by the returned `loglik` closure.
`kerboot()` is unchanged — it resamples rows, and `L`/`D` live outside the data frame.

### Forwarding into primarycensored

In `kerlikelihood()` for `nc == 4`, replace the custom `integrate()` path with a vectorised call to `primarycensored::dprimarycensored()`, grouped by unique `(pwindow, swindow)` pairs.
For each unique `(pwindow_g, swindow_g)` present in `x`, with `pwindow_g = x1r[i] - x1l[i]` and `swindow_g = x2r[i] - x2l[i]`, the group log-contribution is

```
sum( dprimarycensored(
  x         = x2l[group] - x1l[group],
  pdist     = pdist,
  pwindow   = pwindow_g,
  swindow   = swindow_g,
  L         = L,
  D         = D,
  log       = TRUE,
  ...
) )
```

and the total log-likelihood is the sum of group contributions.
`dprimarycensored` dedupes quantile evaluations internally via a lookup table over `sort(unique(c(x, x + swindow)))` (`primarycensored/R/dprimarycensored.R:110-145`), so a single vectorised call is materially faster than the equivalent per-row `pprimarycensored(c(ql, qr), ...)` loop — measured at roughly 5.6x for gamma and 2.7x for gaussian at `n = 100` on integer-day data (see `notes/primarycensored-integration-report.md` section 5).
Grouping by `(pwindow, swindow)` is required because both scalars are positional arguments to `dprimarycensored`; for the common integer-day case with a single primary window of width 1, the group count is 1 and the whole fit is one call.

For `family == "gaussian"`, `pdist = stats::pnorm` with `mean = par1, sd = par2`; for `"gamma"`, `pgamma` with `shape, rate`; for `"lognormal"`, `plnorm`; for `"weibull"`, `pweibull`; for `"skewnorm"`, `pskewnorm` wrapped with `primarycensored::add_name_attribute()` so analytical shortcuts are skipped and numerical integration is used.

`pwindow` must be a positive scalar within each group; rows with `x1l == x1r` (point-observed primary) need a fallback to the single-interval branch below (use the scalar `x1l` as the primary time and evaluate `pdist` directly, with `L`/`D` applied by the same normalisation as 3).

When `L = 0` and `D = Inf`, behaviour must be numerically identical to the current code.
`pprimarycensored` already short-circuits the normalisation step in that case (see `pprimarycensored.R` lines 133–136).
Any remaining difference comes from its integration method vs. the hand-rolled `stats::integrate()` loop; the equivalence test in 5 covers this.

## 2. Semantics

The current per-row log-likelihood is

```
log( F_cens(x2r - x1l) - F_cens(x2l - x1l) )
```

where `F_cens(q) = integral_{0}^{pwindow} F(q - p) f_primary(p) dp` with uniform `f_primary` on `[0, pwindow]`.
This is the probability that a draw from the primary-censored delay distribution lands inside the observed interval `[x2l, x2r]` (re-origined at `x1l`), which is exactly the joint likelihood of the observed primary/secondary windows under uniform primary onset.

Truncation at `[L, D]` is applied to the **underlying delay distribution** (the quantity of interest), not to the observed interval.
With `primarycensored`'s normalisation, `F_cens` is replaced by

```
F_cens_norm(q) = (F_cens(q) - F_cens(L)) / (F_cens(D) - F_cens(L)),   q in [L, D]
F_cens_norm(q) = 0 for q <= L
F_cens_norm(q) = 1 for q >= D
```

and the per-row contribution becomes

```
log( F_cens_norm(x2r - x1l) - F_cens_norm(x2l - x1l) )
  = log( F_cens(x2r - x1l) - F_cens(x2l - x1l) )
  - log( F_cens(D) - F_cens(L) )
```

so each row picks up a `-log(F_cens(D) - F_cens(L))` term.
This is the standard truncated-likelihood correction, summed n times.

The interpretation: data enter the fit as being drawn from the delay distribution conditional on a realised delay falling within `[L, D]`.
The fitted parameters describe the unconditional delay distribution; `[L, D]` is purely an observation filter.

### Boundary subtleties

- `x2r - x1l > D`: the observed interval reaches above the truncation.
  `F_cens_norm(x2r - x1l) = 1` by clamping.
  The row still gets a finite log-lik contribution provided `x2l - x1l < D`.
  We should reject input where `x2l - x1l >= D` because the "real" mass is zero — the row is incompatible with the stated truncation.
- `x2l - x1l < L`: analogous.
  `F_cens_norm(x2l - x1l) = 0` by clamping.
  Reject when `x2r - x1l <= L`.
- `x2l - x1l == L` or `x2r - x1l == D` exactly: see open question below.
  Primarycensored's `.normalise_cdf` takes a `q == L` / `q == D` branch (`pprimarycensored.R` lines 164, 175) that re-uses the unnormalised CDF already computed, so boundary values are handled, but numerical equality may be fragile when `q` comes from floating-point subtraction.

## 3. Single-interval case (nc == 2)

With only `xl, xr` there is no primary window and `primarycensored` is not needed — the underlying delay F itself is what we want.
Proposal: implement truncation inline without calling `pprimarycensored`.
Per row:

```
num  = max(0, pmin(F(xr), F(D)) - pmax(F(xl), F(L)))
den  = F(D) - F(L)
contribution = log(num / den)
```

Equivalent to replacing `F` with the truncated CDF `F_T(q) = (F(min(q,D)) - F(L)) / (F(D) - F(L))` for `q >= L`.
Reject at input check when `xr <= L` or `xl >= D` (row has zero mass under the truncation).

This avoids importing primarycensored's machinery for a single-interval model (no integration is required) and keeps the `nc == 2` path fast.
When `L = 0, D = Inf`, this reduces to the existing `log(F(xr) - F(xl))`.

## 4. Knock-on effects

### `parfitmom` initial values

Moment estimators assume observations come from the untruncated distribution.
Under non-trivial truncation the sample mean and variance of midpoints are biased (shifted toward the centre of `[L, D]`, variance shrunk).
Options:

1. Analytic correction (invert truncated moments). Family-specific, messy, and only yields starting values.
2. Ignore the bias. Starting values are worse but `optim` (Nelder–Mead) usually recovers. Risk: occasional non-convergence on the bootstrap path.
3. Fixed family-dependent starts (e.g. prior `mean ≈ (L + min(D, upperguess)) / 2`).

**Recommendation: option 2.** Accept biased starting points and rely on optim.
Justification: parfitmom is only used to seed `parfitml` and already uses midpoint imputation, so it is already approximate.
Add a comment flagging that under truncation the MoM seed is known-biased.
If bootstrap divergence becomes a problem in testing, add a fallback: try `optim` from the MoM start, then on failure retry from a fixed family-appropriate start before counting a bootstrap discard.

`parfitmom()` still needs to accept `L` and `D` so the signature matches; it can simply ignore them or use them for the fallback start above.

### `kerfeats` mean / variance / quantiles

`parfitml` fits the underlying distribution, so reported features (mean, variance, SD, quantiles) should describe the **untruncated** distribution.
This is what the user wants: `L`, `D` are an observation filter, not part of the estimand.

**Recommendation: leave `kerfeats` reporting untruncated moments and quantiles.**
No code change.
Document in the `parfitml` docstring that features describe the delay distribution as estimated, not the truncated one, and that truncated moments can be recovered externally via `primarycensored::rprimarycensored()` or numerical integration if the user needs them.
Add a note in `summary.parfitml` output when `L > 0` or `D < Inf`, so users see that truncation was applied.

Optional extension (not in this pass): a second feature block in the output called `delayfit_truncated` that holds mean, variance, and quantiles of the truncated distribution, computed from the fitted parameters via `rprimarycensored` sampling or quadrature.

### `kerboot`

No structural change.
`parfitml` already calls `kerboot(x)` then `optim(par, m$loglik, x = xb)`.
Because `L` and `D` are closed over by `m$loglik` at `kerlikelihood` time, they are automatically applied on each bootstrap replicate without `kerboot` knowing about them.

## 5. Test plan

All tests go in a new `tests/testthat/test-truncation.R`.

### Required cases

1. **No-op equivalence**: for each family, simulate double-censored data, fit twice — once with `parfitml(..., L = 0, D = Inf)` (default) and once with the unmodified (pre-integration) fit.
   Assert (a) identical log-likelihood at the MLE to within `1e-6` and (b) identical parameter estimates to within `1e-5`.
   This protects against accidental semantic drift when routing through `primarycensored`.
   Implemented by direct comparison of `m$loglik(v, x)` values at several `v`, not by refitting, so numerical-integration noise is the only source of disagreement.

2. **Parameter recovery under truncation**: simulate `n = 2000` delays from a known gamma (shape = 3, rate = 1), apply uniform primary windows of width 1, discretise the secondary event, filter to `delay in [L = 1, D = 10]`, then fit with `parfitml(..., family = "gamma", L = 1, D = 10)`.
   Assert (a) estimates lie within 10% of truth; (b) unadjusted fit (`L = 0, D = Inf`) on the same truncated sample is materially biased (e.g. mean estimate >15% off truth or outside `parfitml`'s bootstrap 95% CI).
   Run with a fixed `set.seed` for CI determinism.

### Additional cases

3. **Single-interval branch (nc == 2)**: same recovery test with `xl, xr` data (no primary window). Gamma family, `L = 1`, `D = 10`.
4. **Input rejection**: supply a row with `x2r - x1l <= L` or `x2l - x1l >= D` and assert `parfitml` errors with a clear message.
5. **Boundary parity**: a single row where `x2r - x1l == D` to the floating-point digit; assert no `NaN` / `Inf` in the log-lik.
6. **Gaussian with negative `L`**: currently disallowed (`L >= 0`); assert error.
7. **All families smoke test**: fit with `L = 0.5, D = 20` on small simulated data, assert convergence and finite SEs.

## 6. Open questions

- **Boundary equality in `.normalise_cdf`**: `primarycensored/R/pprimarycensored.R:164` uses `any(q == D)` and `any(q == L)`. These are exact float comparisons. When `q = x2r - x1l` and `D` is an integer, subtraction can miss equality. Practical impact is small — it forces an extra `pcens_cdf()` call rather than reusing a cached value — but it may cause tiny log-lik drift. Worth verifying the first no-op test passes to spec; if it doesn't, we may need to round `L`/`D` or document the tolerance.
- **Interval straddling `L` or `D`**: a row with `x2l - x1l < L < x2r - x1l` is legitimate — the observed interval straddles the lower truncation. `pprimarycensored` returns `F_cens_norm(x2l - x1l) = 0` (clamped) and `F_cens_norm(x2r - x1l) = positive`. The difference is the probability mass of the interval above `L`, renormalised. This is correct under the "observation filter on the underlying delay" interpretation, but it does mean a row can contribute non-zero likelihood even though part of its nominal interval is excluded. Confirm this matches the upstream (Gressani) intent before shipping.
- **Row-varying `L`/`D`**: `fitdistdoublecens` accepts per-row `L` and `D`. Real data (e.g. different right-censoring by report date) often need this. Out of scope for the first cut; leave as a follow-up issue.
- **Analytical `pcens_cdf` methods**: `primarycensored` has analytical solutions for some family / primary combinations (see `methods(pcens_cdf)`). For those combinations, the log-lik computation avoids `integrate()` and will be faster than the current code. The no-op test in 5.1 will fail the identity tolerance in those cases — we should relax to `~1e-4` rather than `1e-6`, or force numerical integration for the test.
- **`skewnorm` with primarycensored**: `pskewnorm` is an EpiDelays internal, unknown to `primarycensored`. It must be passed via `add_name_attribute()` so dispatch falls to the numeric path. Verify there is no silent fallback that calls it with unexpected argument names.
- **`parfitmom` under heavy truncation**: if the "accept bias" choice causes observable bootstrap instability in testing, reconsider the fallback-start option.
- **`summary.parfitml` and `kerfeats`**: should summary output surface `L`, `D` explicitly (as fitted arguments) alongside `n` and `Bboot`? Recommended yes, trivial change, not blocking.
- **Vector `pwindow` / `swindow` in `dprimarycensored`**: the design currently groups rows by unique `(pwindow, swindow)` before calling `dprimarycensored`, because both are scalar arguments. Upstream support for vector `pwindow` and `swindow` would eliminate the group-by dance in `kerlikelihood` and let clients pass a single flat vector regardless of window variation. Worth raising as a `primarycensored` feature request; see `notes/primarycensored-integration-report.md` section 6.
