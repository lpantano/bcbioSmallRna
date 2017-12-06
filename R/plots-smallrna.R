#' Small RNA-seq quality control plots
#'
#' @rdname plots-smallrna
#' @author Lorena Pantano, Michael Steinbaugh
#' @aliases bcbSmallMicro bcbSmallCluster bcbSmallSizeDist
#' @param bcb [bcbioSmallRnaDataSet].
#' @param color Column in metadata to use to color the bars.
#' @return [ggplot].
#'
#' @description Plot size distribution of small RNA-seq data.
#' @export
bcbSmallSize <- function(bcb, color = NULL) {
    if (is.null(color))
        color <- metadata(bcb)[["interesting_groups"]]
    info <- adapter(bcb)[["reads_by_sample"]] %>%
        left_join(metrics(bcb), by = "sample")
    m <-  metrics(bcb)
    max_size <- m[["sequence_length"]] %>%
        gsub(".*-", "", .) %>%
        as.numeric() %>%
        .[] / 10L
   m[["library_size"]] <- round(max_size) * 10L

    ggdraw() +
        draw_plot(
            ggplot(info,
                   aes_string(x = color, y = "total", fill = color)) +
                geom_boxplot() +
                ggtitle("total number of reads with adapter") +
                ylab("# reads") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0.5, 1L, 0.5) +
        draw_plot(
            ggplot(m,
                   aes_string(x = color, y = "library_size")) +
                geom_boxplot() +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)),
            0L, 0L, 1L, 0.5)
}

#' @rdname plots-smallrna
#' @keywords plot
#' @export
bcbSmallSizeDist <- function(bcb, color = "NULL"){
    if (is.null(color))
        color <- metadata(bcb)[["interesting_groups"]]

    size <- adapter(bcb)[["reads_by_pos"]] %>%
        left_join(metrics(bcb), by = "sample") %>%
        left_join(adapter(bcb)[["reads_by_sample"]][, c("sample", "total")],
                  by = "sample")
    size[["pct"]] <- size[["V2"]] / size[["total"]] * 100L
    size[["group"]] <- size[[color]]
    ggdraw() +
        draw_plot(
            ggplot(size,
                   aes_string(x = "V1", y = "pct", group = "sample")) +
                geom_smooth(se = FALSE, color = "black") +
                facet_wrap(~group, ncol = 2L) +
                ggtitle("size distribution") +
                ylab("# reads") + xlab("size") +
                theme(
                    axis.text.x = element_text(
                        angle = 90L, vjust = 0.5, hjust = 1L)))
}

#' @rdname plots-smallrna
#' @keywords plot
#' @export
bcbSmallMicro <- function(bcb, group = NULL) {
    if (is.null(group))
        group <- metadata(bcb)[["interesting_groups"]]

    total <- data.frame(sample = colnames(mirna(bcb)),
                        mirtotal = colSums(mirna(bcb)),
                        stringsAsFactors = FALSE) %>%
        left_join(metrics(bcb), by = "sample") %>%
        left_join(adapter(bcb)[["reads_by_sample"]], by ="sample") %>%
        mutate(pct = mirtotal / total * 100L)
    cs <- apply(mirna(bcb), 2L, function(x) {
        cumsum(sort(x, decreasing = TRUE))
    }) %>% as.data.frame %>%
        mutate(pos = 1:nrow(.)) %>%
        melt(., id.vars = "pos") %>%
        left_join(metrics(bcb), by = c("variable" = "sample"))
    es <- melt(mirna(bcb)) %>%
         left_join(metrics(bcb), by = c("X2" = "sample"))
    es[["value"]] <- log2(es[["value"]] + 1L)
    plot_grid(
        ggplot(total) +
            geom_jitter(aes_string(x = group,
                          y = "pct")) +
            ggtitle("Percentage reads being miRNAs") +
            theme(axis.text.x = element_text(
                angle = 90L, hjust = 1L, vjust = 0.5)),
        plot_grid(
            ggplot(es) +
                geom_density(aes_string(x = "value", color = group, group = "X2")) +
                xlab("") +
                ylab("expression") +
                ggtitle("Expression distribution of miRNAs") +
                theme(legend.position="bottom") +
                theme(axis.text.x = element_text(
                    angle = 90L, hjust = 1L, vjust = 0.5)),
            ggplot(cs) +
                geom_line(aes_string(x = "pos",
                                     y = "value",
                                     group = "variable",
                                     color = group)) +
                xlim(0L, 25L) +
                scale_y_log10() +
                theme(legend.position="bottom") +
                ggtitle("Saturation coverage"),
            ncol = 2L),
        nrow = 2L)
}


#' @rdname plots-smallrna
#' @keywords plot
#' @export
bcbSmallCluster <- function(bcb){
    if (is.null(experiments(bcb)[["cluster"]]))
        stop("No cluster data in this analysis.")
    total <- data.frame(sample = colnames(cluster(bcb)),
                        total = colSums(cluster(bcb)))
    ann <- experiments(bcb)[["cluster"]] %>% rowData %>%
        .[["ann"]]
    class <- rep("Other", length(ann))
    class[ann == "|"] <- "None"
    class[grepl("snoRNA", ann, ignore.case = TRUE)] <- "snoRNA"
    class[grepl("repeat", ann, ignore.case = TRUE)] <- "repeat"
    class[grepl("tRNA", ann, ignore.case = TRUE)] <- "tRNA"
    class[grepl("miRNA", ann, ignore.case = TRUE)] <- "miRNA"
    classDF <- bind_cols(bind_cols(cluster(bcb)), ann = class) %>%
        melt() %>% group_by(variable, ann) %>%
        summarise(abundance = sum(value)) %>%
        left_join(., group_by(., variable) %>%
                      summarise(total = sum(abundance))) %>%
        mutate(pct = abundance / total * 100)

    plot_grid(
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
                  axis.text.x = element_text(angle = 90,
                                             vjust = 0.5,
                                             hjust=1)),
        melt(cluster(bcb)) %>%
            ggplot() +
            geom_boxplot(aes_string(x = "variable", y = "value")) +
            xlab("") +
            ylab("expression") +
            scale_y_log10() +
            ggtitle("Expression distribution of clusters detected") +
            theme(axis.text.x = element_text(
            angle = 90L, hjust = 1L, vjust = 0.5)), ncol = 2),
        classDF %>% ggplot(aes_string(x = "variable", y = "pct", fill = "ann")) +
            geom_bar(stat = "identity", position = "dodge") +
            scale_fill_brewer(palette = "Set1") +
            xlab("type of Small RNA") +
            ylab("% of Expression"), nrow = 2
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
#' @param columns Columns in colData to add to heatmap.
#' @param data Only return data.
#' @param ... Options pass to [DEGreport::degPCA()].
#' @export
bcbSmallPCA <- function(bcb, type = "mirna",
                        minAverage = 5, columns = NULL,
                        data = FALSE, ...) {
    if (is.null(columns))
        columns <- metadata(bcb)[["interesting_groups"]]
    counts <- experiments(bcb)[[type]] %>%
        assays %>% .[["rlog"]] %>%
        .[rowMeans(.[]) > minAverage, ]
    annotation_col <- experiments(bcb)[[type]] %>%
        colData %>%
        as.data.frame %>%
        .[, columns, drop = FALSE]
    return(list(counts = counts, annotation = annotation_col))

    th <- HeatmapAnnotation(df = annotation_col)
    hplot <- Heatmap(counts,
            top_annotation = th,
            clustering_method_rows = "ward.D",
            clustering_distance_columns = "kendall",
            clustering_method_columns = "ward.D",
            show_row_names = FALSE,
            show_column_names = ncol(counts) < 50)
    p <- degPCA(counts, annotation_col,
                condition = columns[1], ...)
    print(p)
    hplot
}


