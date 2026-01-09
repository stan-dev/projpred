# Projection predictive feature selection

The R package projpred performs the projection predictive variable (or
"feature") selection for various regression models. We recommend to read
the `README` file (available with enhanced formatting
[online](https://mc-stan.org/projpred/)) and the main vignette
(`topic = "projpred"`, but also available
[online](https://mc-stan.org/projpred/articles/projpred.html)) before
continuing here.

## Terminology

Throughout the whole package documentation, we use the term "submodel"
for all kinds of candidate models onto which the reference model is
projected. For custom reference models, the candidate models don't need
to be actual *sub*models of the reference model, but in any case (even
for custom reference models), the candidate models are always actual
*sub*models of the full
[`formula`](https://rdrr.io/r/stats/formula.html) used by the search
procedure. In this regard, it is correct to speak of *sub*models, even
in case of a custom reference model.

The following model type abbreviations will be used at multiple places
throughout the documentation: GLM (generalized linear model), GLMM
(generalized linear multilevel—or "mixed"—model), GAM (generalized
additive model), and GAMM (generalized additive multilevel—or
"mixed"—model). Note that the term "generalized" includes the Gaussian
family as well.

## Draw-wise divergence minimizers

For the projection of the reference model onto a submodel, projpred
currently relies on the following functions as draw-wise divergence
minimizers (in other words, these are the workhorse functions employed
by projpred's internal default `div_minimizer` functions, see
[`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)):

- Submodel without multilevel or additive terms:

  - For the traditional (or latent) projection (or the augmented-data
    projection in case of the
    [`binomial()`](https://rdrr.io/r/stats/family.html) or
    [`brms::bernoulli()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
    family): An internal C++ function which basically serves the same
    purpose as [`lm()`](https://rdrr.io/r/stats/lm.html) for the
    [`gaussian()`](https://rdrr.io/r/stats/family.html) family and
    [`glm()`](https://rdrr.io/r/stats/glm.html) for all other families.
    The returned object inherits from class `subfit`. Possible tuning
    parameters for this internal C++ function are: `regul` (amount of
    ridge regularization; default: `1e-4`), `thresh_conv` (convergence
    threshold; default: `1e-7`), `qa_updates_max` (maximum number of
    quadratic approximation updates; default: `100`, but fixed to `1` in
    case of the Gaussian family with identity link), `ls_iter_max`
    (maximum number of line search iterations; default: `30`, but fixed
    to `1` in case of the Gaussian family with identity link),
    `normalize` (single logical value indicating whether to scale the
    predictors internally with the returned regression coefficient
    estimates being back-adjusted appropriately; default: `TRUE`),
    `beta0_init` (single numeric value giving the starting value for the
    intercept at centered predictors; default: `0`), and `beta_init`
    (numeric vector giving the starting values for the regression
    coefficients; default: vector of `0`s).

  - For the augmented-data projection:
    [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) (the
    returned object inherits from class `polr`) for the
    [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
    family or
    [`rstanarm::stan_polr()`](https://mc-stan.org/rstanarm/reference/stan_polr.html)
    fits,
    [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
    (the returned object inherits from class `multinom`) for the
    [`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
    family.

- Submodel with multilevel but no additive terms:

  - For the traditional (or latent) projection (or the augmented-data
    projection in case of the
    [`binomial()`](https://rdrr.io/r/stats/family.html) or
    [`brms::bernoulli()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
    family): [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)
    (the returned object inherits from class `lmerMod`) for the
    [`gaussian()`](https://rdrr.io/r/stats/family.html) family,
    [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (the
    returned object inherits from class `glmerMod`) for all other
    families.

  - For the augmented-data projection:
    [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) (the
    returned object inherits from class `clmm`) for the
    [`brms::cumulative()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
    family,
    [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
    (the returned object inherits from class `mmblogit`) for the
    [`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
    family.

- Submodel without multilevel but additive terms:
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) (the returned
  object inherits from class `gam`).

- Submodel with multilevel and additive terms:
  [`gamm4::gamm4()`](https://rdrr.io/pkg/gamm4/man/gamm4.html) (within
  projpred, the returned object inherits from class `gamm4`).

## Verbosity, messages, warnings, errors

Global option `projpred.verbose` may be used for specifying the value
passed to argument `verbose` of
[`project()`](https://mc-stan.org/projpred/dev/reference/project.md),
[`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md), and
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).

By default, projpred catches messages and warnings from the draw-wise
divergence minimizers and throws their unique collection after
performing all draw-wise divergence minimizations (i.e., draw-wise
projections). This can be deactivated by setting global option
`projpred.warn_proj_drawwise` to `FALSE`.

Furthermore, by default, projpred checks the convergence of the
draw-wise divergence minimizers and throws a warning if any seem to have
not converged. This warning is thrown after the warning message from
global option `projpred.warn_proj_drawwise` (see above) and can be
deactivated by setting global option `projpred.check_convergence` to
`FALSE`.

## Parallelization

The projection of the reference model onto a submodel can be run in
parallel (across the projected draws). This is powered by the foreach
package. Thus, any parallel (or sequential) backend compatible with
foreach can be used, e.g., the backends from packages doParallel, doMPI,
or doFuture. Using the global option `projpred.parallel_proj_trigger`,
the number of projected draws below which no parallelization is applied
(even if a parallel backend is registered) can be modified. Such a
"trigger" threshold exists because of the computational overhead of a
parallelization which makes the projection parallelization only useful
for a sufficiently large number of projected draws. By default, the
projection parallelization is turned off, which can also be achieved by
supplying `Inf` (or `NULL`) to option `projpred.parallel_proj_trigger`.
Note that we cannot recommend the projection parallelization on Windows
because in our experience, the parallelization overhead is larger there,
causing a parallel run to take longer than a sequential run. Also note
that the projection parallelization works well for submodels which are
GLMs (and hence also for the latent projection if the submodel has no
multilevel or additive predictor terms), but for all other types of
submodels, the fitted submodel objects are quite big, which—when running
in parallel—may lead to excessive memory usage which in turn may crash
the R session (on Unix systems, setting an appropriate memory limit via
[`unix::rlimit_as()`](https://jeroen.r-universe.dev/unix/reference/rlimit.html)
may avoid crashing the whole machine). Thus, we currently cannot
recommend parallelizing projections onto submodels which are GLMs (in
this context, the latent projection onto a submodel without multilevel
and without additive terms may be regarded as a projection onto a
submodel which is a GLM). However, for
[`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md),
there is also a *CV* parallelization (i.e., a parallelization of
projpred's cross-validation) which can be activated via argument
`parallel` (which in turn can be controlled via global option
`projpred.parallel_cv`).

For the CV parallelization, global option `projpred.export_to_workers`
may be set to a character vector of names of objects to export from the
global environment to the parallel workers.

During parallelization (either of the projection or the CV), progression
updates can be received via the progressr package. This only works if
the doFuture backend is used for parallelization, e.g., via
[`doFuture::registerDoFuture()`](https://doFuture.futureverse.org/reference/registerDoFuture.html)
and `future::plan(future::multisession, workers = 4)`. In that case, the
progressr package can be used, e.g., by calling
`progressr::handlers(global = TRUE)` before running the projection or
the CV in parallel. The projpred package also offers the global option
`projpred.use_progressr` for controlling whether to use the progressr
package (`TRUE` or `FALSE`), but since that global option defaults to
`requireNamespace("progressr", quietly = TRUE) && interactive() && identical(foreach::getDoParName(), "doFuture")`,
it usually does not need to be set by the user.

## Multilevel models: "Integrating out" group-level effects

In case of multilevel models, projpred offers two global options for
"integrating out" group-level effects: `projpred.mlvl_pred_new` and
`projpred.mlvl_proj_ref_new`. When setting `projpred.mlvl_pred_new` to
`TRUE` (default is `FALSE`), then at *prediction* time, projpred will
treat group levels existing in the training data as *new* group levels,
implying that their group-level effects are drawn randomly from a
(multivariate) Gaussian distribution. This concerns both, the reference
model and the (i.e., any) submodel. Furthermore, setting
`projpred.mlvl_pred_new` to `TRUE` causes
[`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
and
[`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md)
to omit the projected group-level effects (for the group levels from the
original dataset). When setting `projpred.mlvl_proj_ref_new` to `TRUE`
(default is `FALSE`), then at *projection* time, the reference model's
fitted values (that the submodels fit to) will be computed by treating
the group levels from the original dataset as *new* group levels,
implying that their group-level effects will be drawn randomly from a
(multivariate) Gaussian distribution (as long as the reference model is
a multilevel model, which—for custom reference models—does not need to
be the case). This also affects the latent response values for a latent
projection correspondingly. Setting `projpred.mlvl_pred_new` to `TRUE`
makes sense, e.g., when the prediction task is such that any group level
will be treated as a new one. Typically, setting
`projpred.mlvl_proj_ref_new` to `TRUE` only makes sense when
`projpred.mlvl_pred_new` is already set to `TRUE`. In that case, the
default of `FALSE` for `projpred.mlvl_proj_ref_new` ensures that at
projection time, the submodels fit to the best possible fitted values
from the reference model, and setting `projpred.mlvl_proj_ref_new` to
`TRUE` would make sense if the group-level effects should be integrated
out completely.

## Memory usage

By setting the global option `projpred.run_gc` to `TRUE`, projpred will
call [`gc()`](https://rdrr.io/r/base/gc.html) at some places (e.g.,
after each size that the forward search passes through) to free up some
memory. These [`gc()`](https://rdrr.io/r/base/gc.html) calls are not
always necessary to reduce the peak memory usage, but they add runtime
(hence the default of `FALSE` for that global option).

## Other notes

Global option `projpred.digits` controls arguments `digits` of
[`print.vselsummary()`](https://mc-stan.org/projpred/dev/reference/print.vselsummary.md)
and
[`print.vsel()`](https://mc-stan.org/projpred/dev/reference/print.vsel.md).

There are several global options to control arguments of
[`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md)
and
[`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md)
globally, see section "Usage" of the help pages of these two functions.

Global option `projpred.warn_L1_interactions` may be set to `FALSE` to
deactivate a warning that an L1 search selected an interaction term
before all involved lower-order interaction terms (including main-effect
terms) were selected (in which case the predictor ranking is
automatically modified by projpred so that the lower-order interaction
terms come before this interaction term).

Most examples are not executed when called via
[`example()`](https://rdrr.io/r/utils/example.html). To execute them,
their code has to be copied and pasted manually to the console.

## Functions

- [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md),
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md):

  For setting up an object containing information about the reference
  model, the submodels, and how the projection should be carried out.
  Explicit calls to
  [`init_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  and
  [`get_refmodel()`](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
  are only rarely needed.

- [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md),
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md):

  For running the *search* part and the *evaluation* part for a
  projection predictive variable selection, possibly with
  cross-validation (CV).

- [`summary.vsel()`](https://mc-stan.org/projpred/dev/reference/summary.vsel.md),
  [`print.vsel()`](https://mc-stan.org/projpred/dev/reference/print.vsel.md),
  [`plot.vsel()`](https://mc-stan.org/projpred/dev/reference/plot.vsel.md),
  [`suggest_size.vsel()`](https://mc-stan.org/projpred/dev/reference/suggest_size.md),
  [`ranking()`](https://mc-stan.org/projpred/dev/reference/ranking.md),
  [`cv_proportions()`](https://mc-stan.org/projpred/dev/reference/cv_proportions.md),
  [`plot.cv_proportions()`](https://mc-stan.org/projpred/dev/reference/plot.cv_proportions.md),
  [`performances()`](https://mc-stan.org/projpred/dev/reference/performances.md):

  For post-processing the results from
  [`varsel()`](https://mc-stan.org/projpred/dev/reference/varsel.md) and
  [`cv_varsel()`](https://mc-stan.org/projpred/dev/reference/cv_varsel.md).

- [`project()`](https://mc-stan.org/projpred/dev/reference/project.md):

  For projecting the reference model onto submodel(s). Typically, this
  follows the variable selection, but it can also be applied directly
  (without a variable selection).

- [`as.matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as.matrix.projection.md)
  and
  [`as_draws_matrix.projection()`](https://mc-stan.org/projpred/dev/reference/as_draws_matrix.projection.md):

  For extracting projected parameter draws.

- [`proj_linpred()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md),
  [`proj_predict()`](https://mc-stan.org/projpred/dev/reference/pred-projection.md):

  For making predictions from a submodel (after projecting the reference
  model onto it).

## See also

Useful links:

- <https://mc-stan.org/projpred/>

- <https://discourse.mc-stan.org>

- Report bugs at <https://github.com/stan-dev/projpred/issues/>

## Author

**Maintainer**: Osvaldo Martin <aloctavodia@gmail.com>

Authors:

- Juho Piironen <juho.t.piironen@gmail.com>

- Markus Paasiniemi

- Alejandro Catalina <alecatfel@gmail.com>

- Frank Weber

- Aki Vehtari

Other contributors:

- Jonah Gabry \[contributor\]

- Marco Colombo \[contributor\]

- Paul-Christian Bürkner \[contributor\]

- Hamada S. Badr \[contributor\]

- Brian Sullivan \[contributor\]

- Sölvi Rögnvaldsson \[contributor\]

- The LME4 Authors (see file 'LICENSE' for details) \[copyright holder\]

- Yann McLatchie \[contributor\]

- Juho Timonen \[contributor\]
