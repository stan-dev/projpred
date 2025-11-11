# Plot ranking proportions from fold-wise predictor rankings

Plots the ranking proportions (see
[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md))
from the fold-wise predictor rankings in a cross-validation with
fold-wise searches. This is a visualization of the *transposed* matrix
returned by
[`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md).
The proportions printed as text inside of the colored tiles are rounded
to whole percentage points (the plotted proportions themselves are not
rounded).

## Usage

``` r
# S3 method for class 'cv_proportions'
plot(
  x,
  text_angle = getOption("projpred.plot_cv_proportions_text_angle", NULL),
  ...
)

# S3 method for class 'ranking'
plot(x, ...)
```

## Arguments

- x:

  For `plot.cv_proportions()`: an object of class `cv_proportions`
  (returned by
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md),
  possibly with `cumulate = TRUE`). For `plot.ranking()`: an object of
  class `ranking` (returned by
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md))
  that
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  will be applied to internally before then calling
  `plot.cv_proportions()`.

- text_angle:

  Passed to argument `angle` of
  [`ggplot2::element_text()`](https://ggplot2.tidyverse.org/reference/element.html)
  for the y-axis tick labels. In case of long predictor names,
  `text_angle = 45` might be helpful (for example).

- ...:

  For `plot.ranking()`: arguments passed to
  [`cv_proportions.ranking()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md)
  and `plot.cv_proportions()`. For `plot.cv_proportions()`: currently
  ignored.

## Value

A ggplot2 plotting object (of class `gg` and `ggplot`).

## Author

Idea and original code by Aki Vehtari. Slight modifications of the
original code by Frank Weber, Yann McLatchie, and Sölvi Rögnvaldsson.
Final implementation in projpred by Frank Weber.

## Examples

``` r
# Data:
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

# The `stanreg` fit which will be used as the reference model (with small
# values for `chains` and `iter`, but only for technical reasons in this
# example; this is not recommended in general):
fit <- rstanarm::stan_glm(
  y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
  QR = TRUE, chains = 2, iter = 1000, refresh = 0, seed = 9876
)

# Run cv_varsel() (with L1 search and small values for `K`, `nterms_max`, and
# `nclusters_pred`, but only for the sake of speed in this example; this is
# not recommended in general):
cvvs <- cv_varsel(fit, method = "L1", cv_method = "kfold", K = 2,
                  nterms_max = 3, nclusters_pred = 10, seed = 5555)
#> Fitting model 1 out of 2
#> Fitting model 2 out of 2

# Extract predictor rankings:
rk <- ranking(cvvs)

# Compute ranking proportions:
pr_rk <- cv_proportions(rk)

# Visualize the ranking proportions:
gg_pr_rk <- plot(pr_rk)
print(gg_pr_rk)


# Since the object returned by plot.cv_proportions() is a standard ggplot2
# plotting object, you can modify the plot easily, e.g., to remove the
# legend:
print(gg_pr_rk + ggplot2::theme(legend.position = "none"))
```
