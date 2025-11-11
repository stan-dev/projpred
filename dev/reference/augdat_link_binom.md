# Link function for augmented-data projection with binomial family

This is the function which has to be supplied to
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
argument `augdat_link` in case of the augmented-data projection for the
[`binomial()`](https://rdrr.io/r/stats/family.html) family.

## Usage

``` r
augdat_link_binom(prb_arr, link = "logit")
```

## Arguments

- prb_arr:

  An array as described in section "Augmented-data projection" of
  [`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
  documentation.

- link:

  The same as argument `link` of
  [`binomial()`](https://rdrr.io/r/stats/family.html).

## Value

An array as described in section "Augmented-data projection" of
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)'s
documentation.
