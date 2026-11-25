# Monotonic effects of `mo()` terms in the parametrization of brms

For a projection with `mo()` terms (see section "Monotonic terms" in
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)),
`mo_effects()` converts the projected coefficients of each `mo()` term
into the parametrization of brms: a coefficient `b` and a simplex
`zeta`, so that the effect of level \\j\\ (counted from 0) is \\b D
\sum\_{k \le j} \zeta_k\\ with \\D\\ the number of levels minus 1.

## Usage

``` r
mo_effects(object, ...)
```

## Arguments

- object:

  An object of class `projection` (returned by
  [`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
  possibly as elements of a `list`).

- ...:

  Currently ignored.

## Value

For a `projection` object, a
[`posterior::draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)
with one row per projected draw (weighted by the weights of the
projected draws) and the columns `bsp_<term>` and `simo_<term>1[k]`
(names as in brms, for example `bsp_mox` and `simo_mox1[1]` for
`mo(x)`). For a `list` of `projection` objects, a `list` of such
matrices.

## Details

For the projected coefficients \\\beta_1, \ldots, \beta_D\\ of the
levels of a `mo()` term, \\b = \sum_k \beta_k / D\\ and \\\zeta_k =
\beta_k / \sum_l \beta_l\\. If all coefficients are 0, `zeta` is `NA`.
In the default monotonic projection (global option
`projpred.mo_monotone`), all coefficients of a term have the same sign,
so `zeta` is on the simplex. In an unconstrained projection, elements of
`zeta` can be negative or larger than 1.

As in brms, an interaction of a `mo()` term with a factor has one
monotonic effect per column of the factor (for example `mo(x):fb`), and
an interaction with a numeric variable has one monotonic effect (for
example `mo(x):z`). Terms with only two levels have no simplex and are
left out, and so are terms with more than one `mo()` call. The
augmented-data projection is not supported. The function needs the
posterior package.
