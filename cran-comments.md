## Test environments

* Local:
    + R version 4.5.1 (2025-06-13) on Ubuntu 24.04.3 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
* win-builder:
    + R-devel R Under development (unstable) (2025-12-04 r89100 ucrt)
    + R-release (R version 4.5.2 (2025-10-31 ucrt))
    + R-oldrelease (version 4.4.3 (2025-02-28 ucrt))


## R CMD check results

All checks gave neither ERRORs nor WARNINGs.

### Local check

The local check on a Linux system gave the following NOTEs:
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/
    
As stated in the first NOTE, the 'cmdstanr' package
(<https://mc-stan.org/cmdstanr/>) is not available on CRAN. The 'cmdstanr'
backend can be used in 'projpred''s unit tests, but it is not necessary for
'projpred' to work. As this NOTE says, we have added the repository URL for the
'cmdstanr' package to the 'DESCRIPTION' file (field 'Additional_repositories').

### win-builder checks

#### R-devel

The R-devel check on win-builder gave the following NOTE:
        
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

For an explanation of this NOTE, see the explanations for the local check above.

#### R-release

The R-release check on win-builder gave the following NOTEs:
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

For an explanation of this NOTE, see the explanations for the local check above.

#### R-oldrelease

The R-oldrelease check on win-builder gave the following NOTEs:

    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

    * checking package dependencies ... NOTE
    Packages suggested but not available for checking: 'unix', 'cmdstanr'

For an explanation of the first NOTE, see the explanations for the local check
above. The unavailability of some "suggested" dependencies mentioned in the
second NOTE is not due to 'projpred'.

## Downstream dependencies

There are two downstream dependencies for this package: 'BayesERtools' and
'brms'. Both of these have been checked locally (with the 'projpred' version
submitted here), at their current CRAN versions and at their most recent
development versions.
