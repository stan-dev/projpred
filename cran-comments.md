## Test environments

* Local:
    + R version 4.1.3 (2022-03-10) on Ubuntu 20.04.4 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
    + R version 4.1.3 (2022-03-10) on Windows 10 x64 (build 19044) system
      (platform: x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
    + R-devel (R version 4.2.0 alpha (2022-03-31 r82049 ucrt))
    + R-release (R version 4.1.3 (2022-03-10))
    + R-oldrelease (R version 4.0.5 (2021-03-31))

## R CMD check results

The local checks gave no ERRORs, WARNINGs, or NOTEs.

For all three win-builder checks, we get a NOTE at
`checking CRAN incoming feasibility ...` concerning the maintainer change (see
above).

Furthermore, for the R-devel and the R-release win-builder checks, we get
another NOTE at `checking CRAN incoming feasibility ...` concerning inaccessible
URLs and DOIs. These URLs are created from DOIs and they are in fact accessible
(just like the DOIs themselves), so we think this can be safely ignored.

## Downstream dependencies

There are two downstream dependencies for this package: 'brms' and 'parameters'.
Both have been checked with respect to the 'projpred' features they use.
