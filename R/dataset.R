#' Tiny example of bcbioSmallRna object
#'
#' Raw fastq file were downloaded from ftp.sra.ebi.ac.uk.
#'
#' Data was analyzed with [bcbio_nextgen](https://bcbio-nextgen.readthedocs.io/en/latest/).
#' Config files to reproduce the analysis are at:
#'
#'   https://github.com/bcbio/bcbio_srnaseq_output_example
#'
#' R code to load the data was:
#'
#' ```
#' sbcb = loadSmallRnaRun("~/repos/bcbio_srnaseq_output_example/final",
#' interestingGroups = "country")
#' ```
#' `
#' @author Lorena Pantano
"sbcb"