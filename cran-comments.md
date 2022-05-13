## Test environments

* Local:
    + R version 4.2.0 (2022-04-22) on Ubuntu 22.04 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
    + R version 4.2.0 (2022-04-22 ucrt) on Windows 10 x64 (build 19044) system
      (platform: x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
    + R-devel (R Under development (unstable) (2022-05-12 r82348 ucrt))
    + R-release (R version 4.2.0 (2022-04-22 ucrt))
    + R-oldrelease (R version 4.1.3 (2022-03-10))

## R CMD check results

The local check gave no ERRORs, WARNINGs, or NOTEs.

For all three win-builder checks, we get a NOTE at
`checking CRAN incoming feasibility ...`:

* On R-devel and R-release, this NOTE reads:
    
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
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

* On R-oldrelease, this NOTE reads:
    
    Maintainer: 'Frank Weber <fweber144[ at ]protonmail.com>'
    
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
    
    Found the following (possibly) invalid DOIs:
      DOI: 10.1214/20-EJS1711
        From: DESCRIPTION
              inst/CITATION
        Status: Internal Server Error
        Message: 500

As can be seen above, the URLs reported as inaccessible are created from DOIs.
In fact, these URLs are accessible (just like the DOIs themselves), so we think
this can be safely ignored.

## Downstream dependencies

There are two downstream dependencies for this package: 'brms' and 'parameters'.
Both have been checked with the 'revdepcheck' package
(<https://r-lib.github.io/revdepcheck>).
