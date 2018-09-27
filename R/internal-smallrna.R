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
#' @param max_samples Maximum samples to perform count transformation.
.read_mirna_counts <- function(meta, col_data, min_hits = 500, max_samples = 50) {
    # TODO Better way to handle sample_dirs than by piping in via metadata?
    fns <- file.path(meta[["sample_dirs"]],
                     paste(names(meta[["sample_dirs"]]),
                           "mirbase-ready.counts",
                           sep = "-"))
    names(fns) <- names(meta[["sample_dirs"]])
    message("Reading miRNA count files")
    fns = .check_compress(fns)
    iso <- IsomirDataSeqFromFiles(
        files = fns[rownames(col_data)],
        coldata = col_data,
        minHits = min_hits,
        design = ~1)
}

.check_compress <- function(fns){
    sapply(fns, function(fn){
        if (file.exists(fn))
            return(fn)
        if (file.exists(paste0(fn, ".gz")))
            return(paste0(fn, ".gz"))
        warning("File doesn't exists: ", fn)
    })
}


.choose_mapping <- function(mapping, priority){
    class <- sapply(priority, function(x){
        grep(x, mapping)
    }) %>% which.min() %>% names(.) %>% .[1L]
    if (is.null(class))
        return(NULL)
    class
}

.clean_id <- function(row){
    priority <- row[1]
    names <- row[2]
    rna <- sapply(str_split(names, ";")[[1]], function(x){
        if (grepl(priority, x))
            return(x)
    })  %>%  unique(.) %>%
        paste(., collapse = ";") %>%
        str_replace_all(., "NULL", "")
    if (rna == "")
        return (str_split(names, ";")[[1]] %>%
                    unique() %>%
                    paste(., collapse = ";"))
    rna
}

#' Load cluster data
#'
#' @keywords internal
#' @rdname read_cluster_counts
#' @author Lorena Pantano
#' @noRd
#' @inheritParams read_smallrna_counts
.read_cluster_counts <- function(meta, col_data, max_samples){
    if (!file.exists(file.path(meta[["project_dir"]], "seqcluster")))
        return(NULL)
    priority = c("miRNA", "tRNA", "repeat", "ncrna", "gene", "rRNA", "")

    clus <- file.path(meta[["project_dir"]],
                     "seqcluster",
                     "counts.tsv") %>%
        read.table(., header = TRUE, sep = "\t",
                   row.names = 1L, check.names = FALSE)
    reads_stats <- file.path(meta[["project_dir"]],
                             "seqcluster",
                             "read_stats.tsv") %>%
        read.table(., header = FALSE, sep = "\t")

    size_data <- file.path(meta[["project_dir"]],
                           "seqcluster",
                           "size_counts.tsv") %>%
        read.table(., header = FALSE, sep = "\t")
    size_data[["cluster"]] <- paste0("cluster:", size_data[["V4"]])
    
    clus_ma <- clus[, 3L:ncol(clus)]
    row_data <- clus[,1L:2L]
    row_data[["cluster"]] <- paste0("cluster:", 1L:nrow(clus))
    row_data <- row_data %>%
        separate(!!sym("ann"), sep = "\\|", into = c("biotype", "names"))
    row_data[["priority"]] <- sapply(row_data[["names"]],
                                     .choose_mapping, priority)
    row_data[["rna_id"]] <- apply(row_data[,c("priority", "names")],
                                  1,
                                  .clean_id)
    row_data <- row_data[,c("cluster", "priority", "rna_id",
                            "biotype", "names", "nloci")]
    row.names(clus_ma) <- paste0("cluster:", 1L:nrow(clus))
    clus_ma <- clus_ma[, row.names(col_data)]
    clus_rlog <- .normalize(clus_ma, col_data, max_samples = max_samples)

    size_data <- left_join(row_data, size_data, by = "cluster") %>% 
        group_by(!!!sym("V1"), !!!sym("priority")) %>% 
        summarise(counts = sum(!!!sym("V2"))) %>% 
        ungroup() %>% 
        mutate(pct = counts / sum(counts) * 100L) %>% 
        dplyr::rename(size = "V1") %>% 
        .[,c("size", "priority", "pct")]
    
    
    SummarizedExperiment(assays = SimpleList(
        raw = clus_ma,
        log = clus_rlog),
        colData = col_data,
        rowData = row_data,
        metadata = list(stats = reads_stats, size = size_data))
}


#' Load cluster data
#'
#' @keywords internal
#' @rdname read_cluster_size
#' @author Lorena Pantano
#' @noRd
#' @inheritParams read_smallrna_counts
.read_cluster_size <- function(meta, col_data, max_samples){
    if (!file.exists(file.path(meta[["project_dir"]], "seqcluster")))
        return(NULL)
    priority = c("miRNA", "tRNA", "repeat", "ncrna", "gene", "")
    
    clus <- file.path(meta[["project_dir"]],
                      "seqcluster",
                      "size_counts.tsv") %>%
        read.table(., header = FALSE, sep = "\t")
    row_data[["priority"]] <- sapply(row_data[["biotype"]],
                                     .choose_mapping, priority)
    row_data[["rna_id"]] <- apply(row_data[,c("priority", "names")],
                                  1,
                                  .clean_id)
    row_data <- row_data[,c("cluster", "priority", "rna_id",
                            "biotype", "names", "nloci")]
    row.names(clus_ma) <- paste0("cluster:", 1L:nrow(clus))
    clus_ma <- clus_ma[, row.names(col_data)]
    clus_rlog <- .normalize(clus_ma, col_data, max_samples = max_samples)
    SummarizedExperiment(assays = SimpleList(
        raw = clus_ma,
        log = clus_rlog),
        colData = col_data,
        rowData = row_data,
        metadata = list(stats = reads_stats))
}

.normalize <- function(ma, col_data, max_samples){
    if (ncol(ma) < max_samples){
        clus_rlog <- ma %>%
            DESeqDataSetFromMatrix(., colData = col_data, design = ~1) %>%
            estimateSizeFactors %>%
            varianceStabilizingTransformation %>% assay
    }else{
        clus_rlog <- voom(ma, plot = FALSE)[["E"]]
    }
    return(clus_rlog)
}

#' Read adapter removal statistics
#'
#' @keywords internal
#' @rdname read_adapter
#' @author Lorena Pantano
#' @noRd
#' @param bcb The [bcbioSmallRnaDataSet] object
.read_adapter <- function(bcb) {
    meta <- metadata(bcb)
    coldata <- colData(bcb)
    fns <- file.path(meta[["sample_dirs"]],
                     paste(names(meta[["sample_dirs"]]),
                           "ready.trimming_stats",
                           sep = "-"))
    names(fns) <- names(meta[["sample_dirs"]])
    reads_by_pos <- lapply(rownames(coldata), function(sample) {
        read.table(fns[sample], sep = "", col.names = c("size", "reads") ) %>%
            mutate(
                sample = sample,
                colorby = coldata[sample, meta[["interesting_groups"]][1]])
    }) %>% bind_rows()
    reads_by_sample <- reads_by_pos %>%
        group_by(!!!sym("sample"), !!!sym("colorby")) %>%
        summarise(total = sum(.data[["reads"]]))
    list(reads_by_pos = reads_by_pos, reads_by_sample = reads_by_sample)
}
