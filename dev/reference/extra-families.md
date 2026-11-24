# Extra family objects

Family objects not in the set of default
[`family`](https://rdrr.io/r/stats/family.html) objects.

## Usage

``` r
Student_t(link = "identity", nu = 3)
```

## Arguments

- link:

  Name of the link function. In contrast to the default
  [`family`](https://rdrr.io/r/stats/family.html) objects, this has to
  be a character string here.

- nu:

  Degrees of freedom for the Student-\\t\\ distribution.

## Value

A family object analogous to those described in
[`family`](https://rdrr.io/r/stats/family.html).

## Details

`Student_t()` is deprecated and will be removed in a future release.
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md), and
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md)
do not support it. For a Student-\\t\\ reference model, use the latent
projection (see
[`extend_family()`](https://mc-stan.org/projpred/dev/reference/extend_family.md)
and section ["Example: Student-t
distribution"](https://mc-stan.org/projpred/articles/latent.html#studentex)
of the latent-projection vignette).
