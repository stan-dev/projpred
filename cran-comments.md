## Test environments

* Local:
    + R version 4.2.2 Patched (2022-11-10 r83330) on Ubuntu 22.04.1 LTS system
      (platform: x86_64-pc-linux-gnu (64-bit))
* win-builder:
    + R-devel (R Under development (unstable) (2023-02-09 r83797 ucrt))
    + R-release (R version 4.2.2 (2022-10-31 ucrt))
    + R-oldrelease (R version 4.1.3 (2022-03-10))

## R CMD check results

All checks gave neither ERRORs nor WARNINGs.

The local check gave the following NOTE:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1002/0471200611
        From: inst/doc/latent.html
        Status: 503
        Message: Service Unavailable

Currently, the 'cmdstanr' package (<https://mc-stan.org/cmdstanr/>) is not
available on CRAN yet. The 'cmdstanr' backend can be used in 'projpred''s unit
tests, but it is not necessary for 'projpred' to work. As this NOTE says, we
have added the repository URL for the 'cmdstanr' package to the 'DESCRIPTION'
file (field 'Additional_repositories').

As can be seen above, the URL reported as inaccessible is created from a DOI. In
fact, this URL is accessible (just like the DOI itself), so we think this can be
safely ignored.

All three win-builder checks gave the following NOTE:

    * checking package dependencies ... NOTE
    Package suggested but not available for checking: 'cmdstanr'

(for an explanation, see above) as well as a NOTE at `checking CRAN incoming
feasibility ...`:

* On R-devel, this latter NOTE reads:
    
    ```
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    ```

* On R-release, this latter NOTE reads:
    
    ```
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1002/0471200611
        From: inst/doc/latent.html
        Status: 503
        Message: Service Unavailable
    ```
    
    For the inaccessible URL, see above.

* On R-oldrelease, this latter NOTE reads:
    
    ```
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Possibly mis-spelled words in DESCRIPTION:
      Bürkner (26:61, 30:29)
      Paasiniemi (25:23)
      Piironen (25:13)
      Vehtari (25:39, 27:5, 28:59, 30:42)
    ```
    
    The "Possibly mis-spelled words in DESCRIPTION" are names that are spelled
    correctly.

The R-release win-builder check also gave the following NOTE:

    * checking examples ... [26s] NOTE
    Examples with CPU (user + system) or elapsed time > 10s
                         user system elapsed
    as.matrix.projection 5.21   0.61   10.48

In most other checks (both local and on win-builder), this example did not take
that long, so the increased runtime might simply be due to a temporary heavy
workload on the system that the win-builder check ran on.

## Downstream dependencies

There is one downstream dependency for this package: 'brms'. This downstream
dependency has been checked with the 'revdepcheck' package
(<https://r-lib.github.io/revdepcheck>) and also locally (with the 'projpred'
version submitted here).
