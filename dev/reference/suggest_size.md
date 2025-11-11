# Suggest submodel size

This function can suggest an appropriate submodel size based on a
decision rule described in section "Details" below. Note that this
decision is quite heuristic and should be interpreted with caution. It
is recommended to examine the results via
[`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md),
[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md),
[`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md),
and/or
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
and to make the final decision based on what is most appropriate for the
problem at hand.

## Usage

``` r
suggest_size(object, ...)

# S3 method for class 'vsel'
suggest_size(
  object,
  stat = "elpd",
  pct = 0,
  type = "upper",
  thres_elpd = NA,
  warnings = TRUE,
  ...
)
```

## Arguments

- object:

  An object of class `vsel` (returned by
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) or
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)).

- ...:

  Arguments passed to
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md),
  except for `object`, `stats` (which is set to `stat`), `type`, and
  `deltas` (which is set to `TRUE`). See section "Details" below for
  some important arguments which may be passed here.

- stat:

  Performance statistic (i.e., utility or loss) used for the decision.
  See argument `stats` of
  [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
  and
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
  for possible choices.

- pct:

  A number giving the proportion (*not* percents) of the *relative* null
  model utility one is willing to sacrifice. See section "Details" below
  for more information.

- type:

  Either `"upper"` or `"lower"` determining whether the decision is
  based on the upper or lower uncertainty interval bound, respectively.
  See section "Details" below for more information.

- thres_elpd:

  Only relevant if `stat %in% c("elpd", "mlpd", "gmpd"))`. The threshold
  for the ELPD difference (taking the submodel's ELPD minus the baseline
  model's ELPD) above which the submodel's ELPD is considered to be
  close enough to the baseline model's ELPD. An equivalent rule is
  applied in case of the MLPD and the GMPD. See section "Details" for a
  formalization. Supplying `NA` deactivates this.

- warnings:

  Mainly for internal use. A single logical value indicating whether to
  throw warnings if automatic suggestion fails. Usually there is no
  reason to set this to `FALSE`.

## Value

A single numeric value, giving the suggested submodel size (or `NA` if
the suggestion failed).

The intercept is not counted by `suggest_size()`, so a suggested size of
zero stands for the intercept-only model.

## Details

In general (beware of special cases below), the suggested model size is
the smallest model size \\j \in \\0, 1, ..., \texttt{nterms\\max}\\\\
for which either the lower or upper bound (depending on argument `type`)
of the uncertainty interval (with nominal coverage `1 - alpha`; see
argument `alpha` of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md))
for \\U_j - U\_{\mathrm{base}}\\ (with \\U_j\\ denoting the \\j\\-th
submodel's true utility and \\U\_{\mathrm{base}}\\ denoting the baseline
model's true utility) falls above (or is equal to) \$\$\texttt{pct}
\cdot (u_0 - u\_{\mathrm{base}})\$\$ where \\u_0\\ denotes the null
model's estimated utility and \\u\_{\mathrm{base}}\\ the baseline
model's estimated utility. The baseline model is either the reference
model or the best submodel found (see argument `baseline` of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)).

In doing so, loss statistics like the root mean squared error (RMSE) and
the mean squared error (MSE) are converted to utilities by multiplying
them by `-1`, so a call such as
`suggest_size(object, stat = "rmse", type = "upper")` finds the smallest
model size whose upper uncertainty interval bound for the *negative*
RMSE or MSE exceeds (or is equal to) the cutoff (or, equivalently, has
the lower uncertainty interval bound for the RMSE or MSE below—or equal
to—the cutoff). This is done to make the interpretation of argument
`type` the same regardless of argument `stat`.

For the geometric mean predictive density (GMPD), the decision rule
above is applied on [`log()`](https://rdrr.io/r/base/Log.html) scale. In
other words, if the true GMPD is denoted by \\U^\ast_j\\ for the
\\j\\-th submodel and \\U^\ast\_{\mathrm{base}}\\ for the baseline model
(so that \\U_j\\ and \\U\_{\mathrm{base}}\\ from above are given by
\\U_j = \log(U^\ast_j)\\ and \\U\_{\mathrm{base}} =
\log(U^\ast\_{\mathrm{base}})\\), then `suggest_size()` yields the
smallest model size whose lower or upper (depending on argument `type`)
uncertainty interval bound for
\\\frac{U^\ast_j}{U^\ast\_{\mathrm{base}}}\\ exceeds (or is equal to)
\$\$(\frac{u^\ast_0}{u^\ast\_{\mathrm{base}}})^{\texttt{pct}}\$\$ where
\\u^\ast_0\\ denotes the null model's estimated GMPD and
\\u^\ast\_{\mathrm{base}}\\ the baseline model's estimated GMPD.

If `!is.na(thres_elpd)` and `stat = "elpd"`, the decision rule above is
extended: The suggested model size is then the smallest model size \\j\\
fulfilling the rule above *or* \\u_j - u\_{\mathrm{base}} \>
\texttt{thres\\elpd}\\. Correspondingly, in case of `stat = "mlpd"` (and
`!is.na(thres_elpd)`), the suggested model size is the smallest model
size \\j\\ fulfilling the rule above *or* \\u_j - u\_{\mathrm{base}} \>
\frac{\texttt{thres\\elpd}}{N}\\ with \\N\\ denoting the number of
observations. Correspondingly, in case of `stat = "gmpd"` (and
`!is.na(thres_elpd)`), the suggested model size is the smallest model
size \\j\\ fulfilling the rule above *or*
\\\frac{u^\ast_j}{u^\ast\_{\mathrm{base}}} \>
\exp(\frac{\texttt{thres\\elpd}}{N})\\.

For example (disregarding the special extensions in case of
`!is.na(thres_elpd)` with `stat %in% c("elpd", "mlpd", "gmpd")`),
`alpha = 2 * pnorm(-1)`, `pct = 0`, and `type = "upper"` means that we
select the smallest model size for which the upper bound of the
`1 - 2 * pnorm(-1)` (approximately 68.3 %) uncertainty interval for
\\U_j - U\_{\mathrm{base}}\\
(\\\frac{U^\ast_j}{U^\ast\_{\mathrm{base}}}\\ in case of the GMPD)
exceeds (or is equal to) zero (one in case of the GMPD), that is (if
`stat` is a performance statistic for which a normal-approximation
uncertainty interval is used, see argument `stats` of
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
and
[`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)),
for which the submodel's utility estimate is at most one standard error
smaller than the baseline model's utility estimate (with that standard
error referring to the utility *difference*).

Apart from the two
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
arguments mentioned above (`alpha` and `baseline`), `resp_oscale` is
another important
[`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md)
argument that may be passed via `...`.

## Examples

``` r
# Data:
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glm(
  y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
  QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
)

# Run varsel() (here without cross-validation, with L1 search, and with small
# values for `nterms_max` and `nclusters_pred`, but only for the sake of
# speed in this example; this is not recommended in general):
vs <- varsel(fit, method = "L1", nterms_max = 3, nclusters_pred = 10,
             seed = 5555)
print(suggest_size(vs))
#> [1] 3
```
