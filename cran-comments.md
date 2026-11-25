## Test environments

* Local:
    + R version 4.5.1 (2025-06-13) on Ubuntu 24.04.3 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))

* win-builder:
    + R-devel: R Under development (unstable) (2026-11-17 r90926 ucrt)
    + R-release: R version 4.6.2 (2026-10-29 ucrt)
    + R-oldrelease: R version 4.5.3 (2026-03-11 ucrt)

## Test environments

* Local:
    + R version 4.5.1 (2025-06-13) on Ubuntu 24.04.3 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))

* win-builder:
    + R-devel: R Under development (unstable) (2026-11-17 r90926 ucrt)
    + R-release: R version 4.6.2 (2026-10-29 ucrt)
    + R-oldrelease: R version 4.5.3 (2026-03-11 ucrt)

## R CMD check results

This is a new major release (3.0.0).

There were no ERRORs or WARNINGs in any of the checks below.

### Local check

The local check on a Linux system gave the following NOTEs:

* Installed package size:

      installed size is 11.8Mb
      sub-directories of 1Mb or more:
        libs  10.0Mb

  The size is due to the compiled code.

* Three further NOTEs are specific to the local system and did not appear on
  win-builder: "unable to verify current time" (future file timestamps), the
  non-portable compilation flag `-mno-omit-leaf-frame-pointer` (from the system's
  R/compiler configuration), and "Skipping checking math rendering: package 'V8'
  unavailable".

### win-builder checks

The win-builder checks printed the following INFO:

    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://stan-dev.r-universe.dev/

The 'cmdstanr' package (<https://mc-stan.org/cmdstanr/>) is not available on CRAN.
The 'cmdstanr' backend can be used in 'projpred''s unit tests, but it is not necessary for 'projpred' to work. The repository URL for 'cmdstanr' is specified in the 'DESCRIPTION' file (field 'Additional_repositories').

#### R-devel

Status: OK

#### R-release

Status: OK

#### R-oldrelease

Status: OK