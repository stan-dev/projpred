## Test environments

* Local:
    + R version 4.5.1 (2025-06-13) on Ubuntu 24.04.3 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))

* win-builder:
    + R-devel R Under development (unstable) (2026-09-29 r90598 ucrt)
    + R-release (R version 4.5.3 (2026-03-11 ucrt))
    + R-oldrelease (R version 4.5.3 (2026-03-11 ucrt))

## R CMD check results

All checks gave neither ERRORs nor WARNINGs.

### Local check

The local check on a Linux system gave the following NOTE:

    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

The 'cmdstanr' package (<https://mc-stan.org/cmdstanr/>) is not available on CRAN.
The 'cmdstanr' backend can be used in 'projpred''s unit tests, but it is not necessary
for 'projpred' to work. The repository URL for 'cmdstanr' is specified in the 'DESCRIPTION' 
file (field 'Additional_repositories').

### win-builder checks

#### R-devel

The R-devel check on win-builder gave the same NOTE regarding 'cmdstanr' as the
local check.

#### R-release

The R-release check on win-builder gave the same NOTE regarding 'cmdstanr' as the
local check.

#### R-oldrelease

The R-oldrelease check on win-builder gave the following NOTEs:

    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

    * checking package dependencies ... NOTE
    Packages suggested but not available for checking: 'unix', 'cmdstanr'

    Possibly misspelled words in DESCRIPTION:
      Paasiniemi
      Piironen
      Vehtari
      rkner

The 'cmdstanr' NOTE is explained above. The unavailability of the suggested
dependencies is not due to 'projpred'. The possible misspellings are proper names
in the package DESCRIPTION and references and are not spelling errors.

## Downstream dependencies

There are two downstream dependencies for this package: 'BayesERtools' and 'brms'. Both of these have been checked locally (with the 'projpred' version submitted here), at their current CRAN versions and at their most recent development versions.