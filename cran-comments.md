## Test environments

* Local:
    + R version 4.5.1 (2025-06-13) on Ubuntu 25.04 system (platform:
      x86_64-pc-linux-gnu (64-bit))
* win-builder:
    + R-devel (R Under development (unstable) (2025-07-02 r88374 ucrt))
    + R-release (R version 4.5.1 (2025-06-13 ucrt))
    + R-oldrelease (R version 4.4.3 (2025-02-28 ucrt))

## R CMD check results

All checks gave neither ERRORs nor WARNINGs.

### Local check

The local check on a Linux system gave the following NOTEs:
    
    * checking CRAN incoming feasibility ... [10s/60s] NOTE
    Maintainer: ‘Osvaldo Martin <[e-mail address]>’
    
    New maintainer:
      Osvaldo Martin <[e-mail address]>
    Old maintainer(s):
      Frank Weber <[e-mail address]>
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/
    
    Found the following (possibly) invalid URLs:
      URL: https://sites.stat.columbia.edu/gelman/arm/
        From: man/mesquite.Rd
        Status: Error
        Message: libcurl error code 60:
          	SSL certificate problem: unable to get local issuer certificate
          	(Status without verification: OK)
      URL: https://sites.stat.columbia.edu/gelman/arm/examples/mesquite/mesquite.dat
        From: man/mesquite.Rd
        Status: Error
        Message: libcurl error code 60:
          	SSL certificate problem: unable to get local issuer certificate
          	(Status without verification: OK)
    
    * checking compilation flags used ... NOTE
    Compilation used the following non-portable flag(s):
      ‘-mno-omit-leaf-frame-pointer’
    
    * checking examples ... [27s/28s] NOTE
    Examples with CPU (user + system) or elapsed time > 5s
                          user system elapsed
    as.matrix.projection 8.788  0.437   9.376

As stated in the first NOTE, we would like to change the maintainer in this
release.

As also stated in the first NOTE, the 'cmdstanr' package
(<https://mc-stan.org/cmdstanr/>) is not available on CRAN. The 'cmdstanr'
backend can be used in 'projpred''s unit tests, but it is not necessary for
'projpred' to work. As this NOTE says, we have added the repository URL for the
'cmdstanr' package to the 'DESCRIPTION' file (field 'Additional_repositories').

The URLs reported as possibly invalid in the first NOTE are in fact accessible.
The certificate problems are not related to 'projpred'.

The compilation flag reported in the second NOTE is not set within 'projpred',
so it is probably due to the system that this local check ran on. In CRAN
checks, this NOTE has never appeared yet (but it is shown when checking the
previous CRAN version 2.8.0 locally).

The runtime of the 'as.matrix.projection' example (see the third NOTE) may
sometimes be slightly > 5s.

### win-builder checks

#### R-devel

The R-devel check on win-builder gave the following NOTE:
    
    * checking CRAN incoming feasibility ... [15s] NOTE
    Maintainer: 'Osvaldo Martin <[e-mail address]>'
    
    New maintainer:
      Osvaldo Martin <[e-mail address]>
    Old maintainer(s):
      Frank Weber <[e-mail address]>
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

For an explanation of this NOTE, see the explanations for the local check above.

#### R-release

The R-release check on win-builder gave the following NOTEs:
    
    * checking CRAN incoming feasibility ... [14s] NOTE
    Maintainer: 'Osvaldo Martin <[e-mail address]>'
    
    New maintainer:
      Osvaldo Martin <[e-mail address]>
    Old maintainer(s):
      Frank Weber <[e-mail address]>
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

For an explanation of this NOTE, see the explanations for the local check above.

#### R-oldrelease

The R-oldrelease check on win-builder gave the following NOTEs:
    
    * checking CRAN incoming feasibility ... [20s] NOTE
    Maintainer: 'Osvaldo Martin <[e-mail address]>'
    
    New maintainer:
      Osvaldo Martin <[e-mail address]>
    Old maintainer(s):
      Frank Weber <[e-mail address]>
    
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
