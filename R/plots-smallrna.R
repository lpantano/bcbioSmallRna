#' Small RNA-seq quality control plots
#'
#' @rdname plots-smallrna
#' @author Lorena Pantano, Michael Steinbaugh
#' @aliases bcbSmallMicro bcbSmallCluster
#' @param bcb [bcbioSmallRnaDataSet].
#' @param color Column in metadata to use to color the bars.
#' @return [ggplot].
#'
#' @description Plot size distribution of small RNA-seq data.
#' @export
bcbSmallSize <- function(bcb, color="group") {
    info <- adapter(bcb)
    ggdraw() +
        draw_plot(
            ggplot(info[["reads_by_sample"]],
                   aes_string(x = "sample", y = "total", fill = color)) +
                geom_bar(stat = "identity", position = "dodge") +
                ggtitle("total number of reads with adapter") +
                ylab("# reads") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0.5, 1L, 0.5) +
        draw_plot(
            ggplot(info[["reads_by_pos"]],
                   aes_string(x = "V1", y = "V2", group = "sample")) +
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
#' @keywords plot
#' @export
bcbSmallMicro <- function(bcb) {
    total <- data.frame(sample = colnames(mirna(bcb)),
                     total = colSums(mirna(bcb)))
    cs <- apply(mirna(bcb), 2L, function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }) %>% as.data.frame
    cs[["pos"]] <- 1L:nrow(cs)
    plot_grid(
        ggplot(total) +
            geom_bar(aes_string(x = "sample",
                          y = "total"),
                     stat = "identity") +
            ggtitle("Total reads being miRNAs") +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(mirna(bcb))) +
            geom_boxplot(aes_string(x = "X2", y = "value")) +
            xlab("") +
            ylab("expression") +
            scale_y_log10() +
            ggtitle("Expression distribution fo miRNAs") +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        ggplot(melt(cs, id.vars = "pos")) +
            geom_line(aes_string(x = "pos",
                           y = "value",
                           color = "variable")) +
            xlim(0L, 50L) +
            scale_y_log10() +
            ggtitle("Saturation coverage"),
        nrow = 3L)
}


#' @rdname plots-smallrna
#' @keywords plot
#' @export
bcbSmallCluster <- function(bcb){
    if (is.null(experiments(bcb)[["cluster"]]))
        stop("No cluster data in this analysis.")
    total <- data.frame(sample = colnames(cluster(bcb)),
                        total = colSums(cluster(bcb)))
    plot_grid(
        experiments(bcb)[["cluster"]] %>%
            metadata %>%
            .[["stats"]] %>%
            ggplot(aes_string(x = "V2", y = "V1", fill = "V3")) +
            geom_bar(stat = 'identity', position = 'dodge') +
            labs(list(x="samples", y="reads")) +
            scale_fill_brewer("steps", palette = 'Set1') +
            ggtitle("Reads kept after filtering") +
            theme(legend.position = "bottom",
                  legend.direction = "horizontal",
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
        melt(cluster(bcb)) %>%
            ggplot() +
            geom_boxplot(aes_string(x = "variable", y = "value")) +
            xlab("") +
            ylab("expression") +
            scale_y_log10() +
            ggtitle("Expression distribution of clusters detected") +
            theme(axis.text.x = element_text(
            angle = 90L, hjust = 1L, vjust = 0.5))
    )
}

#' Plot PCA of the different small RNA data
#' @rdname plotsSmallPCA
#' @author Lorena Pantano
#' @keywords plot
#' @description Clustering small RNA samples.
#' @param bcb [bcbioSmallRnaDataSet].
#' @param type Data type to plot: `mirna, cluster, isomir, trna`.
#' @param minAverage Minimun average small RNA expression to be kept.
#' @param ... Options pass to [DEGreport::degPCA()].
#' @export
bcbSmallPCA <- function(bcb, type = "mirna", minAverage = 5, ...) {
    counts <- experiments(bcb)[[type]] %>%
        assays %>% .[["rlog"]] %>%
        .[colMeans(.[]) > 2, ]
    annotation_col <- colData(bcb) %>%
        as.data.frame %>%
        .[, metadata(bcb)[["interesting_groups"]], drop = FALSE]
    th <- HeatmapAnnotation(df = annotation_col)
    hplot <- Heatmap(counts,
            top_annotation = th,
            clustering_method_rows = "ward.D",
            clustering_distance_columns = "kendall",
            clustering_method_columns = "ward.D",
            show_row_names = FALSE)
    p <- degPCA(counts, annotation_col,
                condition = metadata(bcb)[["interesting_groups"]], ...)
    print(p)
    hplot
}


