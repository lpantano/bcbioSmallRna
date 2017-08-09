#' Create isomiRs object from bcbio output
#'
#' Read bcbio sample information from YAML to get isomiR object.
#'
#' @rdname read_smallrna_counts
#' @keywords internal
#'
#' @author Lorena Patano
#'
#' @param meta metadata of the bcbio run.
.read_smallrna_counts <- function(meta) {
    # TODO Better way to handle sample_dirs than by piping in via metadata?
    fns <- file.path(meta[["sample_dirs"]],
                     paste(names(meta[["sample_dirs"]]),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(meta[["sample_dirs"]])
    message("Reading miRNA count files")
    iso <- IsomirDataSeqFromFiles(
        files = fns[rownames(meta[["metadata"]])],
        coldata = meta[["metadata"]],
        design = ~meta[["interesting_groups"]][1L])
    isoNorm(iso)
}
