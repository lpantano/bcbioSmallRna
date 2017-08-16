#' Create isomiRs object from bcbio output
#'
#' Read bcbio sample information from YAML to get isomiR object.
#'
#' @rdname read_smallrna_counts
#' @keywords internal
#'
#' @author Lorena Pantano
#' @noRd
#' @param meta Metadata of the bcbio run.
#' @param col_data Samples information.
.read_mirna_counts <- function(meta, col_data) {
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

#' Load cluster data
#'
#' @keywords internal
#' @rdname read_cluster_counts
#' @author Lorena Pantano
#' @noRd
#' @inheritParams read_smallrna_counts
.read_cluster_counts <- function(meta, col_data){
    if (!file.exists(file.path(meta[["project_dir"]], "seqcluster")))
        return(NULL)
    clus <- file.path(meta[["project_dir"]],
                     "seqcluster",
                     "counts.tsv") %>%
        read.table(., header = TRUE, sep = "\t", row.names = 1L)
    reads_stats <- file.path(meta[["project_dir"]],
                             "seqcluster",
                             "read_stats.tsv") %>%
        read.table(., header = FALSE, sep = "\t")

    ann <- clus[, 2L]
    clus_ma <- clus[, 3L:ncol(clus)]
    row.names(clus_ma) <- paste0(row.names(clus_ma), ann)
    clus_ma <- clus_ma[, row.names(col_data)]
    clus_rlog <- clus_ma %>%
        DESeqDataSetFromMatrix(., colData = col_data, design = ~1) %>%
        estimateSizeFactors %>%
        rlog %>% assay
    SummarizedExperiment(assays = SimpleList(
        raw = clus_ma,
        rlog = clus_rlog),
        colData = col_data,
        rowData = clus[,1L:2L],
        metadata = list(stats = reads_stats))
}

#' Read adapter removal statistics
#'
#' @keywords internal
#' @rdname read_adapter
#' @author Lorena Pantano
#' @noRd
#' @param bcb [bcbioSmallRnaDataSet]
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
        group_by(!!!sym("sample"), !!!sym("group")) %>%
        summarise(total = sum(.data[["V2"]]))
    list(reads_by_pos = reads_by_pos, reads_by_sample = reads_by_sample)
}
