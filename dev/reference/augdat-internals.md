# Augmented-data projection: Internals

The augmented-data projection makes extensive use of *augmented-rows
matrices* and *augmented-length vectors*. In the following, \\N\\,
\\C\_{\mathrm{cat}}\\, \\C\_{\mathrm{lat}}\\, \\S\_{\mathrm{ref}}\\, and
\\S\_{\mathrm{prj}}\\ from help topic
[refmodel-init-get](https://mc-stan.org/projpred/dev/reference/refmodel-init-get.md)
are used. Furthermore, let \\C\\ denote either \\C\_{\mathrm{cat}}\\ or
\\C\_{\mathrm{lat}}\\, whichever is appropriate in the context where it
is used (e.g., for `ref_predfun`'s output, \\C = C\_{\mathrm{lat}}\\).
Similarly, let \\S\\ denote either \\S\_{\mathrm{ref}}\\ or
\\S\_{\mathrm{prj}}\\, whichever is appropriate in the context where it
is used. Then an augmented-rows matrix is a matrix with \\N \cdot C\\
rows in \\C\\ blocks of \\N\\ rows, i.e., with the \\N\\ observations
nested in the \\C\\ (possibly latent) response categories. For ordered
response categories, the \\C\\ (possibly latent) response categories
(i.e., the row blocks) have to be sorted increasingly. The columns of an
augmented-rows matrix have to correspond to the \\S\\ parameter draws,
just like for the traditional projection. An augmented-rows matrix is of
class `augmat` (inheriting from classes `matrix` and `array`) and needs
to have the value of \\C\\ stored in an attribute called `ndiscrete`. An
augmented-length vector (class `augvec`) is the vector resulting from
subsetting an augmented-rows matrix to extract a single column and
thereby dropping dimensions.
