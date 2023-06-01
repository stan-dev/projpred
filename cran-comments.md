## Test environments

* Local:
    + R version 4.2.3 (2023-03-15) on Ubuntu 22.04.2 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
    + R version 4.3.0 (2023-04-21) on Ubuntu 22.04.2 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
* win-builder:
    + R-devel (R Under development (unstable) (2023-05-30 r84466 ucrt))
    + R-release (R version 4.3.0 (2023-04-21 ucrt))
    + R-oldrelease (R version 4.2.3 (2023-03-15 ucrt))
* macOS builder:
    + R version 4.3.0 Patched (2023-05-18 r84451) on macOS Ventura 13.3.1 system
      (platform: aarch64-apple-darwin20 (64-bit))
      (r-release-macosx-arm64|4.3.0|macosx|macOS 13.3.1 (22E261)|Mac mini|Apple
      M1||en_US.UTF-8|macOS 11.3|clang-1403.0.22.14.1|GNU Fortran (GCC) 12.2.0)

## R CMD check results

All checks gave neither ERRORs nor WARNINGs.

The local check with R version 4.3.0 gave the following NOTE:
    
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/

The local check with R version 4.2.3 gave the following NOTE:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1002/0471200611
        From: inst/doc/latent.html
        Status: 403
        Message: Forbidden

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
    ```

* On R-oldrelease, this latter NOTE reads:
    
    ```
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1002/0471200611
        From: inst/doc/latent.html
        Status: 403
        Message: Forbidden
    ```
    
    For the inaccessible URL, see above.

The macOS builder check gave the following NOTE:

    * checking package dependencies ... NOTE
    Package suggested but not available for checking: 'cmdstanr'

(for an explanation, see above).

## Downstream dependencies

There is one downstream dependency for this package: 'brms'. This downstream
dependency has been checked with the 'revdepcheck' package
(<https://r-lib.github.io/revdepcheck>) and also locally (with the 'projpred'
version submitted here).
