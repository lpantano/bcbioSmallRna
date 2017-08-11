#' Small RNA-seq quality control plots
#'
#' @rdname plots-smallrna
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @param bcb [bcbioSmallRnaDataSet].
#' @return [ggplot].
#'
#' @description Plot size distribution of small RNA-seq data.
#' @export
plot_size_distribution <- function(bcb) {
    info <- adapter(bcb)
    ggdraw() +
        draw_plot(
            ggplot(info[["reads_by_sample"]],
                   aes_(x = ~sample, y = ~total, fill = ~group)) +
                geom_bar(stat = "identity", position = "dodge") +
                ggtitle("total number of reads with adapter") +
                ylab("# reads") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0.5, 1L, 0.5) +
        draw_plot(
            ggplot(info[["reads_by_pos"]], aes_(x = ~V1, y = ~V2, group = ~sample)) +
                geom_bar(stat = "identity", position = "dodge") +
                facet_wrap(~group, ncol = 2L) +
                ggtitle("size distribution") +
                ylab("# reads") + xlab("size") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0L, 1L, 0.5)
}



#' @rdname plots-smallrna
#' @description Plot of total miRNA counts.
#' @export
plot_mirna_counts <- function(bcb) {
    .t <- data.frame(sample = colnames(bcbio(bcb)),
                     total = colSums(bcbio(bcb)))
    cs <- apply(bcbio(bcb), 2L, function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }) %>% as.data.frame
    cs[["pos"]] <- 1L:nrow(cs)
    plot_grid(
        ggplot(.t) +
            geom_bar(aes_(x = ~sample,
                          y = ~total),
                     stat = "identity") +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(bcbio(bcb))) +
            geom_boxplot(aes_(x = ~X2, y = ~value)) +
            xlab("") +
            ylab("expression") +
            scale_y_log10() +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(cs, id.vars = "pos")) +
            geom_line(aes_(x = ~pos,
                           y = ~value,
                           color = ~variable)) +
            xlim(0L, 50L) +
            scale_y_log10(),
        nrow = 3L)
}


#' @rdname plots-smallrna
#' @description Clustering small RNA samples.
#' @export
plot_srna_clusters <- function(bcb) {
    counts <- bcbio(bcb)
    design <- metadata(bcb)[["metadata"]]
    dds <- DESeqDataSetFromMatrix(
        counts[rowSums(counts > 0L) > 3L, ],
        colData = design,
        design = ~1L)
    vst <- rlog(dds, betaPriorVar = FALSE)
    annotation_col <- design %>%
        .[, metadata(bcb)[["interesting_groups"]], drop = FALSE]
    pheatmap(assay(vst),
             annotation_col = annotation_col,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D",
             scale = "row",
             show_rownames = FALSE)
    plot_pca(metadata(bcb), vst)
}
