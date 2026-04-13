# Design: Non-uniform primary event distributions via primarycensored

Scope: add a user-facing primary event density argument to `kerlikelihood()` and `parfitml()` so that the doubly interval-censored engine can fit data generated under non-uniform primary onset (e.g. exponential growth).
Upstream support already lives in `primarycensored::dprimarycensored(..., dprimary, dprimary_args)`; the EpiDelays work is forwarding, validation, and surfacing.
This note is design-only; no R source, tests, or other non-`notes/` files are changed.
The sister design note `notes/truncation-design.md` covers `L`/`D` truncation; this one assumes that feature is being added in the same series and specifies how the two compose.

## 1. API proposal

### Naming: `dprimary` vs `primary_dist`

The EpiDelays public surface uses short, unpunctuated names (`family`, `Bboot`, `pgbar`, and in the truncation design `L`, `D`).
The existing `build_pc_loglik()` helper in `R/kerlikelihood.R` already binds directly to primarycensored's argument layout via `pdist`, `pwindow`, `swindow`, `log`.
Calling the new argument `dprimary` mirrors primarycensored exactly, keeps the forwarding logic trivial, and surfaces the same vocabulary users will encounter when they read `?primarycensored::dprimarycensored`.
It also avoids a second spelling that implementers have to keep in sync if primarycensored evolves.

**Recommendation: `dprimary`, companion args `dprimary_args`.**
This applies to both `kerlikelihood()` and `parfitml()`.

### Companion args format

`primarycensored::dprimarycensored` takes `dprimary_args = list(...)` and passes the list through to the primary density via `do.call()`.
The three candidate styles on the EpiDelays side are:

1. `dprimary_args = list(r = 0.1)` — mirrors upstream exactly.
2. `...` captured on `parfitml()` and re-split into delay vs primary args.
3. A named list inlined into a wrapper closure the user builds (`dprimary = function(x, ...) dexpgrowth(x, r = 0.1, ...)`).

Option 3 works today without any API change but shifts boilerplate onto the user and silently breaks analytical-solution dispatch unless the user remembers `add_name_attribute()`.
Option 2 collides with the existing `...` surface and makes validation messages ambiguous.
Option 1 is the most forgiving, matches upstream, and lets EpiDelays stay a thin forwarder.

**Recommendation: `dprimary_args = list()`, forwarded verbatim.**

### String shorthand

Accepting `dprimary = "uniform"` or `dprimary = "expgrowth"` is friendlier but adds a lookup table that has to be maintained against `primarycensored::pcd_primary_distributions` and obscures the namespace of the actual function being called.
Raw functions are also what primarycensored itself exposes.
A beginner who already knows `stats::dunif` can type it directly; a beginner who doesn't will usually land on the docstring, which can point at `primarycensored::dexpgrowth`.

**Recommendation: accept a function only; no string shorthand in this pass.**
If demand appears later, a thin wrapper can be added without breaking the function-accepting surface.

### Signatures

```r
kerlikelihood(
  x,
  family,
  likapprox     = "ni",
  L             = 0,
  D             = Inf,
  dprimary      = stats::dunif,
  dprimary_args = list()
)

parfitml(
  x,
  family,
  Bboot         = 1000,
  pgbar         = TRUE,
  L             = 0,
  D             = Inf,
  dprimary      = stats::dunif,
  dprimary_args = list()
)
```

`parfitmom()` does not gain `dprimary` (see 6).

### Roxygen `@param` entries

```
#' @param dprimary Primary event density function. Used only when \code{x} has
#' four columns (doubly interval-censored data). Must take a vector \code{x}
#' and the arguments \code{min} and \code{max}, return a density normalised to
#' integrate to 1 on \code{[min, max]}, and be compatible with
#' \code{primarycensored::dprimarycensored()}. Defaults to \code{stats::dunif}
#' (uniform primary onset within the observation window), which reproduces the
#' behaviour of earlier EpiDelays versions. Non-uniform choices include
#' \code{primarycensored::dexpgrowth} for exponential growth during an
#' outbreak. Ignored when \code{x} has two columns (single interval-censored
#' data).
#'
#' @param dprimary_args A named list of additional arguments passed to
#' \code{dprimary}, mirroring the \code{dprimary_args} argument of
#' \code{primarycensored::dprimarycensored()}. Example:
#' \code{list(r = 0.1)} for \code{primarycensored::dexpgrowth}. Defaults to
#' an empty list. Must be empty unless \code{x} has four columns.
```

## 2. Math

With a non-uniform primary density `f_primary` on `[0, pwindow]`, the primary-censored CDF becomes

```
F_cens(q; pwindow) = integral_0^pwindow F(q - p) f_primary(p) dp
```

and the per-row log-likelihood contribution is

```
log( F_cens(x2r - x1l; pwindow) - F_cens(x2l - x1l; pwindow) )
```

exactly as in the uniform case but without the `1/pwindow` factor — that factor was an artefact of `f_primary` being uniform, and it drops out when `f_primary` is supplied directly.
`primarycensored::dprimarycensored` already performs this integral internally once `dprimary` and `dprimary_args` are passed in, so the EpiDelays side is a one-line forwarding change inside `build_pc_loglik()`: the two new arguments flow through the `do.call()` alongside the existing `pdist`, `pwindow`, `swindow`, `log` slots.
No further per-row arithmetic is required.

## 3. Single-interval branch (`nc == 2`)

Single interval-censored data have no primary event window; `dprimary` is meaningless there.
The cleanest contract is to accept the argument uniformly (so a caller can thread a single `parfitml(...)` call through both branches without conditional argument assembly) and validate that it has not been overridden when `nc == 2`.

**Recommendation: validate at entry.** In `kerlikelihood()` and `parfitml()`, when `nc == 2`:

1. If `dprimary` is not the default `stats::dunif` (identity-checked with `identical(dprimary, stats::dunif)`), error with a clear message: "`dprimary` only applies to doubly interval-censored data (four columns); `x` has two columns so no primary event window is modelled".
2. If `dprimary_args` is non-empty (`length(dprimary_args) > 0`), error with the same message.
3. Otherwise ignore both and proceed on the existing `nc == 2` code path.

Identity-check against `stats::dunif` rather than checking the default at the call site is slightly looser but survives users writing `dprimary = stats::dunif` explicitly, which is fine.

## 4. Interaction with the mc approximation

`kerlikelihood(..., likapprox = "mc")` currently draws primary onsets via `stats::runif(n = M, min = x$x1l[i], max = x$x1r[i])` inside each family's nc == 4 loop (`R/kerlikelihood.R:139`, `:179`, `:219`, `:257`, `:295`).
A non-uniform `dprimary` makes the uniform sampler wrong: the Monte Carlo estimate of `F_cens` would be biased, silently.

Three options:

(a) Support non-uniform primary on `likapprox = "ni"` only; error cleanly if the user combines `likapprox = "mc"` with a non-default `dprimary`.
(b) Support both by asking the user for a companion `rprimary` sampler (primarycensored ships `rexpgrowth`) and swapping the `runif()` call for it in each family branch.
(c) Defer the mc branch entirely as out of scope.

Option (b) doubles the argument surface (`dprimary`, `dprimary_args`, `rprimary`, `rprimary_args`) for a code path that `?kerlikelihood` already describes as slow and fallback-only.
Option (c) leaves an unused fallback alive without touching it, but implicitly re-defines what the mc branch means (uniform-primary fallback even when the ni branch has been told otherwise) and will surprise users.

**Recommendation: option (a).**
On the mc branch, if `identical(dprimary, stats::dunif)` and `length(dprimary_args) == 0` the current code runs unchanged; otherwise error with "`likapprox = \"mc\"` does not support non-uniform `dprimary`; use `likapprox = \"ni\"` (the default)".
This is cheap, preserves the mc branch for existing users, and does not silently produce biased estimates.

## 5. Interaction with truncation

The truncation design adds `L` and `D` to `kerlikelihood()` and `parfitml()` and routes them into `primarycensored::dprimarycensored(..., L, D, ...)`.
Composing with `dprimary` is order-preserving: primarycensored first builds `F_cens` using the supplied `dprimary` (section 2 above), then applies `.normalise_cdf(L, D)` on top to produce `F_cens_norm`.
EpiDelays therefore just forwards both sets of arguments into the same `do.call()` and the engine composes them correctly.

The `build_pc_loglik()` call inside `kerlikelihood()` becomes, for a single `(pwindow, swindow)` group:

```r
logd <- do.call(
  primarycensored::dprimarycensored,
  c(
    list(
      x             = lowers[idx],
      pdist         = pdist,
      pwindow       = pw,
      swindow       = sw,
      L             = L,
      D             = D,
      dprimary      = dprimary,
      dprimary_args = dprimary_args,
      log           = TRUE
    ),
    pars
  )
)
```

All of `L`, `D`, `dprimary`, `dprimary_args` are closed over at `kerlikelihood()` time and live on the `loglik` closure, so `parfitml()` and its bootstrap loop do not need to know about them individually.

Row-level validation ordering (before the first `loglik` evaluation):

1. `primarycensored::.check_truncation_bounds(L, D)`.
2. `primarycensored::check_dprimary(dprimary, pwindow, dprimary_args)` per unique `pwindow` in the data (`check_dprimary` integrates the primary density over `[0, pwindow]` to confirm it is normalised).
3. EpiDelays row-range checks against `L`/`D` as per `notes/truncation-design.md` section 1.

`check_dprimary` additionally enforces that `dprimary` takes `min` and `max` arguments, which rules out user closures that hardcode a pwindow and rules in `stats::dunif` / `primarycensored::dexpgrowth`.
This matches the upstream contract, so any function that works with `dprimarycensored` will work here.

## 6. Knock-on effects

### `parfitmom`

Moment matching assumes observations come from an untruncated, uniform-primary model.
Under a non-uniform primary the midpoint sample mean is shifted (`dexpgrowth` with `r > 0` weights late-window primary onsets, so midpoints overshoot the delay mean by a family-dependent amount) and the variance is deflated slightly.
Analytic correction is family- and primary-distribution-specific and adds complexity out of proportion to the benefit.

**Recommendation: leave `parfitmom` unchanged and unaware of `dprimary`.**
Justification mirrors the truncation design note section 4.1: MoM values seed `optim` and Nelder-Mead recovers from moderate bias.
`parfitml()` does not forward `dprimary` to `parfitmom()`, and no new argument is added to `parfitmom()`.
Add a comment in `parfitml.R` flagging that the MoM seed is known-biased under non-uniform primary, and (if bootstrap discards become noisy in the test suite) fall back to a fixed family-appropriate start on bootstrap optim failure.

### `kerfeats`

Features describe the underlying delay distribution parameters, not the primary window.
No change needed.

### `kerboot`

Unchanged.
The bootstrap resamples rows of `x`; `dprimary` and `dprimary_args` live on the closure `m$loglik` and are automatically applied on every replicate.

### `summary.parfitml`

`summary.parfitml` (`R/S3_summary_parfitml.R`) prints a block of fixed rows (family, npars, censtype, n, elapsed, mleconv, Bboot, discards, AIC, BIC) in both `full` and `compact` modes.
Add a new row: `Primary event dist :` rendered as the function's name (or "uniform (default)" when `identical(object$dprimary, stats::dunif)`).
If `dprimary_args` is non-empty, append a one-line summary of the list entries, e.g. `r = 0.1`.
This keeps the surface consistent with the truncation design, which surfaces `L`, `D`, `n`, and `Bboot` on the same block.
Storing the two new entries on the `parfitml` object is a matter of propagating them from `kerlikelihood()`'s return list through `parfitml()`'s `o <- c(...)` assembly.

## 7. Test plan

Tests go in a new `tests/testthat/test-dprimary.R` (paralleling `test-truncation.R`).

1. **No-op equivalence.** For each family (`gaussian`, `gamma`, `lognormal`, `weibull`, `skewnorm`), call
   `parfitml(x, family, dprimary = stats::dunif, dprimary_args = list())`
   and assert (a) identical numeric log-likelihood at several parameter points to the default-constructed fit (which does not pass `dprimary`), and (b) identical MLE parameter estimates to within `1e-6`.
   Protects against forwarding bugs.
2. **Non-uniform primary recovery.** Simulate doubly-censored gamma data with known parameters using `primarycensored::rprimarycensored(..., dprimary = primarycensored::dexpgrowth, dprimary_args = list(r = 0.1))`, discretise to unit-width primary and secondary windows, fit with the matching `dprimary = primarycensored::dexpgrowth, dprimary_args = list(r = 0.1)`, and assert parameter estimates are within 10% of truth with `n = 2000` and a fixed `set.seed`.
3. **Mismatched primary is biased.** Same simulated data as (2), fit with the default uniform, assert that the `shape` or `rate` estimate is outside a 10% window of truth.
   Catches the silent-failure mode where a user forgets to specify `dprimary` on non-uniform data.
4. **Single-interval rejection.** `kerlikelihood(x_nc2, family = "gamma", dprimary = primarycensored::dexpgrowth, dprimary_args = list(r = 0.1))` errors with a message naming `dprimary`.
   Same for `parfitml`.
   Also assert that `dprimary = stats::dunif` and `dprimary_args = list()` do not error on nc == 2 data.
5. **End-to-end with truncation.** `parfitml(x, family = "gamma", L = 1, D = 10, dprimary = primarycensored::dexpgrowth, dprimary_args = list(r = 0.1))` converges on simulated data (and the same `r` used to simulate).
   Covers the composition of sections 4 and 5.
6. **mc branch rejection.** `kerlikelihood(x, family = "gamma", likapprox = "mc", dprimary = primarycensored::dexpgrowth, dprimary_args = list(r = 0.1))` errors clearly.
   `kerlikelihood(x, family = "gamma", likapprox = "mc")` still works (default uniform).
7. **Invalid `dprimary`.** Passing a non-normalised function errors via `primarycensored::check_dprimary`; check that the error is surfaced cleanly and references the user-facing argument name, not an internal one.

## 8. Open questions

- **Cost of custom `dprimary` with no analytical `pcens_cdf` method.** Dispatch falls through to `pcens_cdf.default`, which uses numerical integration per unique quantile. This is the same path that the skew-normal row in `notes/primarycensored-integration-report.md` section 3 table already incurs, so the per-call overhead is a known quantity. For non-uniform primary on gamma/lognormal/weibull, the analytical shortcuts tagged `pcens_pgamma_dunif` / `pcens_plnorm_dunif` / `pcens_pweibull_dunif` stop applying (they are `_dunif`-suffixed), so every family collapses onto the numeric default. Fit times should be benchmarked after implementation; a rough expectation is that gamma with `dexpgrowth` will be roughly `gaussian + dunif` cost, since both lose the analytical path.
- **Re-exporting `primarycensored::dexpgrowth`.** For ergonomics, a user-facing `parfitml(..., dprimary = dexpgrowth, dprimary_args = list(r = 0.1))` is nicer than `primarycensored::dexpgrowth`. A re-export costs one line in `NAMESPACE` and one `#' @importFrom primarycensored dexpgrowth` tag but commits EpiDelays to any upstream rename. **Recommendation: do not re-export in this pass; document the fully qualified form and revisit if users request it.**
- **Vignette / roxygen example.** The introductory `parfitml()` example in roxygen should stay on the default uniform case to keep the primary example minimal. Add a separate `\dontrun{}` block showing `primarycensored::dexpgrowth` with a non-trivial `r`, so the feature is discoverable from the manpage without complicating the first-read example.
- **Interaction with a future `pskewnorm` with non-uniform primary.** `pskewnorm` is passed into `dprimarycensored` via the `pskewnorm_cdf` wrapper in `build_pc_loglik()` without `add_name_attribute()`, so it already routes through `pcens_cdf.default` regardless of `dprimary`. There is nothing to change, but the test matrix in 7.1 should include `skewnorm + dunif` to catch any regression from the forwarding change.
- **Row-varying `dprimary`.** Real data can have time-varying primary onset (e.g. growth rate changing between epochs). Out of scope for this pass; scalar `dprimary` shared across all rows only. Raise as a follow-up issue if a user asks, and consider whether it fits cleanly into the row-varying `L`/`D` extension listed in `notes/truncation-design.md` section 6.
- **LOW CONFIDENCE: `check_dprimary` under multiple unique `pwindow` values.** `primarycensored::check_dprimary(dprimary, pwindow, dprimary_args)` takes a scalar `pwindow` and integrates over `[0, pwindow]`. When `x` contains several distinct primary window widths, EpiDelays should call `check_dprimary` once per unique width before entering the optim loop; the alternative is to let `dprimarycensored` run its own check per group on every likelihood evaluation, which is wasteful. Worth confirming the upstream check is cheap enough that the latter is acceptable, or wrapping it in EpiDelays entry-time validation.
