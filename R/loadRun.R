#' Load small RNA-seq bcbio-nextgen run
#'
#' Simply point to the final upload directory output by
#' [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/), and this function
#' will take care of the rest. It automatically imports small RNA-seq counts,
#' metadata, and program versions used.
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param interestingGroups Character vector of interesting groups. First entry
#'   is used for plot colors during quality control (QC) analysis. Entire vector
#'   is used for PCA and heatmap QC functions.
#' @param maxSamples *Optional*. Maximum number of samples to calculate rlog
#'   and variance stabilization object from DESeq2.
#' @param dataDir Folder to keep a cache of the object.
#' @param ... Additional arguments, saved as metadata.
#'
#' @return [bcbioSmallRnaDataSet].
#' @importFrom yaml yaml.load_file
#' @export
loadSmallRnaRun <- function(
    uploadDir = "final",
    interestingGroups = "description",
    maxSamples = 50,
    dataDir = NULL,
    colData = NULL,
    ...) {
    message("Cache will be safed under ", dataDir)
    # Directory paths and cache path====
    if (!is.null(dataDir)) {
        if (file.exists(file.path(dataDir, "bcb.rda"))){
            load(file.path(dataDir, "bcb.rda"))
            return(bcb)
        }
    }

    if (!dir.exists(uploadDir)) {
        stop("Final upload directory failed to load")
    }
    upload_dir <- normalizePath(uploadDir)
    # Find most recent nested project_dir (normally only 1)
    project_dir_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
    project_dir <- dir(upload_dir,
                       pattern = project_dir_pattern,
                       full.names = FALSE,
                       recursive = FALSE)
    if (length(project_dir) != 1L) {
        stop("Uncertain about project directory location")
    }
    message(project_dir)
    match <- str_match(project_dir, project_dir_pattern)
    run_date <- match[[2L]] %>% as.Date
    template <- match[[3L]]
    project_dir <- file.path(upload_dir, project_dir)


    # Project summary YAML ====
    yaml_file <- file.path(project_dir, "project-summary.yaml")
    if (!file.exists(yaml_file)) {
        stop("YAML project summary missing")
    }
    message("Reading project summary YAML")
    yaml <- yaml.load_file(yaml_file)

    # colData ====
    if (is.null(colData)){
        col_data <- .yaml_metadata(yaml)
    }else{
        col_data <- as.data.frame(colData)
        col_data[["sample"]] <- rownames(col_data)
    }


    # Sample names ====
    # Obtain the samples (and their directories) from the YAML
    sample_names <- vapply(
        yaml[["samples"]],
        function(x) x[["description"]],
        character(1L)) %>% sort %>%
        intersect(., rownames(col_data))
    sample_dirs <- file.path(upload_dir, sample_names) %>%
        set_names(sample_names)
    if (!identical(basename(sample_dirs), sample_names)) {
        stop("Sample name assignment mismatch")
    }
    message(paste(length(sample_dirs), "samples detected"))


    # Genome ====
    # Use the genome build of the first sample to match
    genome_build <- yaml[["samples"]][[1L]][["genome_build"]]
    message(paste("Genome:", genome_build))

    # Sequencing lanes ====
    lane_pattern <- "_L(\\d{3})"
    if (any(str_detect(sample_dirs, lane_pattern))) {
        lanes <- str_match(names(sample_dirs), lane_pattern) %>%
            .[, 2L] %>%
            unique %>%
            length
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"))
    } else {
        # TODO Check that downstream functions don't use NULL
        lanes <- 1L
    }


    # Sample metrics ====
    metrics <- .yaml_metrics(yaml) %>%
        .[.[["description"]] %in% col_data[["sample"]],]

    # bcbio-nextgen run information ====
    message("Reading bcbio run information")
    data_versions <- .data_versions(project_dir)
    programs <- .programs(project_dir)
    bcbio_nextgen <- read_lines(
        file.path(project_dir, "bcbio-nextgen.log"))
    bcbio_nextgen_commands <- read_lines(
        file.path(project_dir, "bcbio-nextgen-commands.log"))


    # Metadata ====
    metadata <- list(
        analysis = "smallRna",
        upload_dir = upload_dir,
        sample_dirs = sample_dirs,
        project_dir = project_dir,
        template = template,
        run_date = run_date,
        interesting_groups = interestingGroups,
        genome_build = genome_build,
        lanes = lanes,
        yaml_file = yaml_file,
        yaml = yaml,
        metrics = metrics,
        data_versions = data_versions,
        programs = programs,
        bcbio_nextgen = bcbio_nextgen,
        bcbio_nextgen_commands = bcbio_nextgen_commands)

    # SummarizedExperiment for miRNA ====
    mirna <- .read_mirna_counts(metadata, col_data, max_samples = maxSamples)
    mirna_rlog <- mirna %>%
        isoNorm(., maxSamples = maxSamples) %>% counts(., norm = TRUE)
    isomirna <- isoCounts(mirna, iso5 = TRUE, iso3 = TRUE,
                        add = TRUE, subs = TRUE, seed = TRUE,
                        ref = TRUE)
    iso_rlog <- isomirna %>%
        isoNorm(., maxSamples = maxSamples) %>% counts(., norm = TRUE)
    mir <- SummarizedExperiment(assays = SimpleList(
        raw = counts(mirna),
        rlog = mirna_rlog),
        colData = col_data[rownames(colData(mirna)),, FALSE])
    iso <- SummarizedExperiment(assays = SimpleList(
        raw = counts(isomirna),
        rlog = iso_rlog),
        colData = col_data[rownames(colData(isomirna)),, FALSE])
    # SummarizedExperiment for tRNA ====

    # SummarizedExperiment for clusters ====
    cluster <- .read_cluster_counts(metadata, col_data, max_samples = maxSamples)
    # SummarizedExperiment for mirdeep2 ====

    # MultiAssayExperiment ====
    exps <- list(mirna = mir, isomirs = iso, cluster = cluster)
    se <- .multiassay_experiment(exps, col_data, metadata)

    # bcbioSmallRnaDataSet ====
    bcb <- new("bcbioSmallRnaDataSet", se)
    bcbio(bcb, "isomirs") <- mirna
    adapter(bcb) <- .read_adapter(bcb)
    if (!is.null(dataDir)){
        if(!file.exists(dataDir))
            dir.create(dataDir, showWarnings = FALSE, recursive = TRUE)
        save(bcb, file = file.path(dataDir, "bcb.rda"))
    }
    bcb
}
