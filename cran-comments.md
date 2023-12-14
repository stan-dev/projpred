## Test environments

* Local:
    + R version 4.3.2 (2023-10-31) on Ubuntu 22.04.3 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
    + R version 4.2.3 (2023-03-15 ucrt) on Windows 10 system (platform:
      x86_64-w64-mingw32 (64-bit))
* win-builder:
    + R-devel (R Under development (unstable) (2023-12-13 r85679 ucrt))
    + R-release (R version 4.3.2 (2023-10-31 ucrt))
    + R-oldrelease (R version 4.2.3 (2023-03-15 ucrt))
* macOS builder:
    + R version 4.3.0 Patched (2023-05-18 r84451) on macOS Ventura 13.3.1 system
      (platform: aarch64-apple-darwin20 (64-bit))
      (r-release-macosx-arm64|4.3.0|macosx|macOS 13.3.1 (22E261)|Mac mini|Apple
      M1||en_US.UTF-8|macOS 11.3|clang-1403.0.22.14.1|GNU Fortran (GCC) 12.2.0)

## R CMD check results

Except for the R-oldrelease check on win-builder (see below), all checks gave
neither ERRORs nor WARNINGs.

### Local checks

#### Linux

The local check with R version 4.3.2 under Linux gave the following NOTE:
    
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/

Currently, the 'cmdstanr' package (<https://mc-stan.org/cmdstanr/>) is not
available on CRAN yet. The 'cmdstanr' backend can be used in 'projpred''s unit
tests, but it is not necessary for 'projpred' to work. As this NOTE says, we
have added the repository URL for the 'cmdstanr' package to the 'DESCRIPTION'
file (field 'Additional_repositories').

#### Windows

The local check with R version 4.2.3 under Windows gave the following NOTEs:
    
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
    
    * checking package dependencies ... NOTE
    Package suggested but not available for checking: 'unix'
    
    * checking Rd cross-references ... NOTE
    Package unavailable to check Rd xrefs: 'unix'
    
    * checking examples ... [18s] NOTE
    Examples with CPU (user + system) or elapsed time > 5s
                         user system elapsed
    as.matrix.projection 5.45    0.3    5.75

As can be seen in the first NOTE, the URL reported as inaccessible is created
from a DOI. In fact, this URL is accessible, so we think this can be safely
ignored. For an explanation of the 'cmdstanr'-related content of the first NOTE,
see above.

The 'unix' package mentioned in the second and the third NOTE is not available
for Windows, but it is only mentioned in 'projpred''s documentation, so has no
impact on functionality.

The runtime of the 'as.matrix.projection' example (see the fourth NOTE) may
sometimes be slightly > 10s (see also below).

### win-builder checks

#### R-devel

The R-devel check on win-builder gave the following NOTEs:
    
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    * checking package dependencies ... NOTE
    Package suggested but not available for checking: 'unix'

For an explanation of these NOTEs, see the explanations for the local checks
above.

#### R-release

The R-release check on win-builder gave the following NOTEs:
    
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    * checking package dependencies ... NOTE
    Packages suggested but not available for checking: 'unix', 'cmdstanr'
    
    * checking examples ... [38s] NOTE
    Examples with CPU (user + system) or elapsed time > 10s
                          user system elapsed
    as.matrix.projection 11.02   0.53   11.52

For an explanation of these NOTEs, see the explanations for the local checks
above.

#### R-oldrelease

The R-oldrelease check on win-builder gave the following NOTEs:
    
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
    
    * checking package dependencies ... NOTE
    Packages suggested but not available for checking: 'unix', 'cmdstanr'
    
    * checking examples ... [42s] NOTE
    Examples with CPU (user + system) or elapsed time > 10s
                          user system elapsed
    as.matrix.projection 12.21   0.53   12.76

For an explanation of these NOTEs, see the explanations for the local checks
above.

The R-oldrelease check on win-builder also gave the following ERROR:
    
    Running the tests in 'tests/testthat.R' failed.
    [...]
    Caused by error in `initializePtr()`:
      ! function 'cholmod_factor_ldetA' not provided by package 'Matrix'
    [...]
    lme4 (local) initializePtr()

However, this error was not reproducible on a local Windows system with R 4.2.3
(see the local checks above). Furthermore, this error does not seem to be
related to 'projpred', but rather at least one of its dependencies. Hence, there
is probably nothing that can be done in 'projpred' to avoid this error.

### macOS builder

The macOS builder check gave the following NOTE:
    
    * checking package dependencies ... NOTE
    Package suggested but not available for checking: 'cmdstanr'

For an explanation of this NOTE, see the explanations for the local checks
above.

## Downstream dependencies

There is one downstream dependency for this package: 'brms'. This downstream
dependency has been checked locally (with the 'projpred' version submitted
here).
