## Test environments

* Local:
    + R version 4.2.1 (2022-06-23) on Ubuntu 22.04.1 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
    + R version 4.2.1 (2022-06-23 ucrt) on Windows 10 x64 (build 19044) system
      (platform: x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
    + R-devel (R Under development (unstable) (2022-08-17 r82724 ucrt))
    + R-release (R version 4.2.1 (2022-06-23 ucrt))
    + R-oldrelease (R version 4.1.3 (2022-03-10))
* macOS builder:
    + R version 4.2.1 Patched (2022-06-23 r82516) on macOS 11.5.2 (20G95) system
      (platform: aarch64-apple-darwin20 (64-bit))
      (r-release-macosx-arm64|4.2.1|macosx|macOS 11.5.2 (20G95)|Mac mini|Apple
      M1||en_US.UTF-8)

## R CMD check results

All checks gave neither ERRORs nor WARNINGs.

The local checks gave the following NOTE:

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

All three win-builder checks gave the following NOTE:

    * checking package dependencies ... NOTE
    Package suggested but not available for checking: 'cmdstanr'

(for an explanation, see above) as well as a NOTE at `checking CRAN incoming
feasibility ...`:

* On R-devel, this latter NOTE reads:
    
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1214/17-EJS1337SI
        From: inst/doc/projpred.html
        Status: 500
        Message: Internal Server Error
      URL: https://doi.org/10.1214/20-EJS1711
        From: inst/doc/projpred.html
              README.md
        Status: 500
        Message: Internal Server Error

* On R-release, this latter NOTE reads:
    
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Suggests or Enhances not in mainstream repositories:
      cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1214/17-EJS1337SI
        From: inst/doc/projpred.html
        Status: 500
        Message: Internal Server Error
      URL: https://doi.org/10.1214/20-EJS1711
        From: inst/doc/projpred.html
              README.md
        Status: 500
        Message: Internal Server Error
    
    Found the following (possibly) invalid DOIs:
      DOI: 10.1214/20-EJS1711
        From: DESCRIPTION
              inst/CITATION
        Status: Internal Server Error
        Message: 500

* On R-oldrelease, this latter NOTE reads:
    
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
    Possibly mis-spelled words in DESCRIPTION:
      Bürkner (22:41)
      Paasiniemi (21:27)
      Piironen (21:17)
      Vehtari (21:42, 22:53)
    
    Found the following (possibly) invalid URLs:
      URL: https://doi.org/10.1214/17-EJS1337SI
        From: inst/doc/projpred.html
        Status: 500
        Message: Internal Server Error
      URL: https://doi.org/10.1214/20-EJS1711
        From: man/projpred-package.Rd
              inst/doc/projpred.html
              README.md
        Status: 500
        Message: Internal Server Error

As can be seen above, the URLs reported as inaccessible are created from DOIs.
In fact, these URLs are accessible (just like the DOIs themselves), so we think
this can be safely ignored.

The "Possibly mis-spelled words in DESCRIPTION" are names that are spelled
correctly.

The macOS builder check gave the following NOTE:

    * checking package dependencies ... NOTE
    Package suggested but not available for checking: ‘cmdstanr’

## Downstream dependencies

There is one downstream dependency for this package: 'brms'. This downstream
dependency has been checked with the 'revdepcheck' package
(<https://r-lib.github.io/revdepcheck>).
