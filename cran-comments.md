## Maintainer change

With this submission, we would like to change the 'projpred' maintainer from
Alejandro Catalina (<alecatfel@gmail.com>) to Frank Weber
(<fweber144@protonmail.com>). Alejandro has sent an e-mail to
<CRAN-submissions@R-project.org> on March 31, 2022, stating that he agrees with
the maintainer change.

## Test environments

* Local:
    + R 4.1.3 on Ubuntu 20.04.3 LTS system (platform:
      x86_64-pc-linux-gnu (64-bit))
    + R 4.1.3 on Windows 10 x64 (build 19044) system (platform:
      x86_64-w64-mingw32/x64 (64-bit))
* win-builder:
    + R-devel (R version 4.2.0 alpha (2022-03-24 r81978 ucrt))
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

Additionally, the R-devel win-builder check gave an ERROR because of a unit test
failure. In local tests, we have never gotten this failure before, so we hope it
will not occur anymore in future R-devel versions or in the final R 4.2.0
version.

## Downstream dependencies

There are two downstream dependencies for this package: 'brms' and 'parameters'.
Both have been checked with respect to the 'projpred' features they use. While
doing so, we discovered some long-standing 'projpred'-related bugs in the
'parameters' package which we have reported on the 'parameters' issue tracker.
