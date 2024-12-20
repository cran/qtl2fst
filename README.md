# qtl2fst

[![R-CMD-check](https://github.com/rqtl/qtl2fst/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rqtl/qtl2fst/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/qtl2fst)](https://cran.r-project.org/package=qtl2fst)
[![zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3237776.svg)](https://doi.org/10.5281/zenodo.3237776)

[R/qtl2](https://kbroman.org/qtl2/) (aka qtl2) is a reimplementation of
the QTL analysis software [R/qtl](https://rqtl.org), to better handle
high-dimensional data and complex cross designs.

The [qtl2fst](https://github.com/rqtl/qtl2fst) package uses
the [fst package](https://www.fstpackage.org/) to store genotype probabilities in
a set of files for rapid access but reduced memory usage. It is fully
integrated with [R/qtl2](https://kbroman.org/qtl2/), and is modeled after the
[qtl2feather package](https://github.com/byandell-sysgen/qtl2feather) which stores
files in the [feather](https://github.com/wesm/feather) format.

---

### Installation

Install R/qtl2fst from [CRAN](https://cran.r-project.org):

    install.packages("qtl2fst")

---

### Vignette

View the [qtl2fst user
guide](https://kbroman.org/qtl2/assets/vignettes/qtl2fst.html),
available [online](https://kbroman.org/qtl2/assets/vignettes/qtl2fst.html)
and also within R:

    vignette("qtl2fst", package="qtl2fst")



---

#### License

Licensed under [GPL-3](https://www.r-project.org/Licenses/GPL-3).
