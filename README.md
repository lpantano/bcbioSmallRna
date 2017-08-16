# bcbioSmallRna


[![Build Status](https://travis-ci.org/lpantano/bcbioSmallRna.svg?branch=master)](https://travis-ci.org/lpantano/bcbioSmallRna)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)


Quality control and differential expression for [bcbio][] small RNA-seq experiments.


## Installation

This is an [R][] package.

### [Bioconductor][] method

```r
source("https://bioconductor.org/biocLite.R")
biocLite("lpantano/bcbioSmallRna")
```

### [devtools][] method

```r
install.packages("devtools")
devtools::install_github("lpantano/bcbioSmallRna")
```


[bcbio]: https://github.com/chapmanb/bcbio-nextgen
[bioconductor]: https://bioconductor.org
[devtools]: https://cran.r-project.org/package=devtools
[r]: https://www.r-project.org
[rmarkdown]: http://rmarkdown.rstudio.com