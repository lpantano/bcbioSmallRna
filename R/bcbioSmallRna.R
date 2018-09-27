#' bcbioSmallRna
#'
#' Quality control and differential expression for bcbio-nextgen small RNA-seq
#' experiments.
#'
#' @import BiocGenerics Biobase DESeq2 MultiAssayExperiment SummarizedExperiment
#'   S4Vectors isomiRs
#' @importFrom limma voom
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom DEGreport degQC degCovariates degPatterns degPCA
#' @importFrom ggplot2 aes_string coord_fixed coord_flip coord_trans element_blank element_text
#'   expand_limits facet_wrap geom_bar geom_boxplot geom_density geom_hline
#'   geom_jitter geom_line geom_point geom_polygon geom_ribbon geom_smooth
#'   ggplot ggtitle guides labs scale_x_continuous scale_y_log10 theme xlab xlim
#'   ylab ylim scale_x_log10 scale_color_manual geom_text aes_string
#'   scale_x_log10 annotation_logticks scale_fill_brewer
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr select mutate mutate_if filter left_join right_join
#'                   arrange "%>%"
#'                   bind_rows bind_cols mutate_all summarise group_by
#'                   enquo inner_join ungroup
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom janitor remove_empty
#' @importFrom reshape melt
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr separate
#' @importFrom purrr set_names
#' @importFrom isomiRs IsomirDataSeqFromFiles
#' @importFrom methods as is new slot slot<- validObject
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom stats formula
#' @importFrom stringr str_detect str_match str_split str_replace_all
#' @importFrom readr read_lines read_csv read_tsv read_delim
#' @importFrom rlang sym dots_values
#' @importFrom utils read.table capture.output download.file
"_PACKAGE"

globalVariables(".")

# Quality control plot colors
fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"

# Plot label separator
label_sep <- " : "
