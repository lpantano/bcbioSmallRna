#' Download Dependency File
#'
#' If the required dependency file isn't present, download latest version from
#' the [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' File download utility for RMarkdown knit reports.
#'
#' @rdname download
#' @name download
#'
#' @examples
#' bcbSmallDownload()
NULL



# Constructors ====
.download <- function(file) {
    sapply(seq_along(file), function(a) {
        if (!file.exists(file[[a]])) {
            download.file(
                file.path("https://raw.githubusercontent.com/lpantano/bcbioSmallRna/master/docs",
                          "downloads",
                          file[[a]]),
                destfile = file[[a]])
        }
    }) %>% invisible
}



# Methods ====
#' @rdname download
#' @export
bcbSmallDownload <- function() {
    .download(
        c("_output.yaml",
          "_footer.Rmd",
          "_header.Rmd",
          "setup.R"))
}
