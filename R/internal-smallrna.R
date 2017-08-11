#' Create isomiRs object from bcbio output
#'
#' Read bcbio sample information from YAML to get isomiR object.
#'
#' @rdname read_smallrna_counts
#' @keywords internal
#'
#' @author Lorena Patano
#' @noRd
#' @param meta metadata of the bcbio run.
.read_smallrna_counts <- function(meta, col_data) {
    # TODO Better way to handle sample_dirs than by piping in via metadata?
    fns <- file.path(meta[["sample_dirs"]],
                     paste(names(meta[["sample_dirs"]]),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(meta[["sample_dirs"]])
    message("Reading miRNA count files")
    iso <- IsomirDataSeqFromFiles(
        files = fns[rownames(col_data)],
        coldata = col_data,
        design = ~1)
    isoNorm(iso)
}

.read_adapter <- function(bcb) {
    meta <- metadata(bcb)
    coldata <- colData(bcb)
    fns <- file.path(meta[["sample_dirs"]],
                     paste(names(meta[["sample_dirs"]]),
                           "ready.trimming_stats",
                           sep = "-"))
    names(fns) <- names(meta[["sample_dirs"]])
    reads_by_pos <- lapply(rownames(coldata), function(sample) {
        read.table(fns[sample], sep = "") %>%
            mutate(
                sample = sample,
                group = coldata[sample, meta[["interesting_groups"]]])
    }) %>% bind_rows()
    reads_by_sample <- reads_by_pos %>%
        group_by(!!!sym("sample") , !!!sym("group")) %>%
        summarise(total = sum(.data[["V2"]]))
    list(reads_by_pos = reads_by_pos, reads_by_sample = reads_by_sample)
}
