.make_names <- function(names) {
    make.names(names, unique = TRUE) %>%
        gsub("[^[:alnum:] ]", "_", .) %>%
        gsub("\\.", "_", .) %>%
        tolower
}

#' YAML utilities
#'
#' @rdname yaml
#' @keywords internal
#' @author Michael Steinbaugh, Lorena Pantano
#'
#' @param yaml Project summary YAML.
#' @param keys Nested operator keys supplied as dot objects.
#'
#' @note Metrics are only generated for a standard RNA-seq run with aligned
#'   counts. Fast RNA-seq mode with lightweight counts (pseudocounts) doesn't
#'   output the same metrics into the YAML.
#' @noRd
#' @importFrom rlang dots_values
#' @return [DataFrame].
.yaml <- function(yaml, keys) {
    samples <- yaml[["samples"]]
    if (!length(samples)) {
        stop("No sample information in YAML")
    }

    # Check for nested keys, otherwise return NULL
    # TODO Improve the recursion method using sapply in a future update
    if (!keys[[1L]] %in% names(samples[[1L]])) {
        return(NULL)
    }
    if (length(keys) > 1L) {
        if (!keys[[2L]] %in% names(samples[[1L]][[keys[[1L]]]])) {
            return(NULL)
        }
    }

    m <- lapply(seq_along(samples), function(a) {
        nested <- samples[[a]][[keys]] %>%
            set_names(., nm = .make_names(names(.)))
        # Set the description
        nested[["description"]] <- samples[[a]][["description"]]
        # Remove legacy duplicate `name` identifier
        nested[["name"]] <- NULL
        if (rev(keys)[[1L]] == "metadata") {
            if (is.null(nested[["batch"]]))
                nested[["batch"]] <- NA
            if (is.null(nested[["phenotype"]]))
                nested[["phenotype"]] <- NA
            if (length(nested[["phenotype"]])) {
                if (grepl("^$", nested[["phenotype"]]))
                    nested[["phenotype"]] <- NA
            }
        }
        nested[!sapply(nested, is.null)] %>%
            data.frame(., stringsAsFactors = FALSE)
    }
    ) %>%
        bind_rows %>%
        arrange(description) %>%
        as.data.frame %>%
        remove_empty_cols
    m
}



#' @rdname yaml
.yaml_metadata <- function(yaml) {
    metadata <- .yaml(yaml, "metadata")
    if (ncol(metadata) == 1)
        metadata[["group"]] <- metadata[["description"]]
    metadata %>%
        as.data.frame %>%
        mutate_all(factor) %>%
        mutate(description = as.character(.data[["description"]])) %>%
        as.data.frame %>%
        column_to_rownames("description") %>%
        DataFrame
}



#' @rdname yaml
.yaml_metrics <- function(yaml) {
    metrics <- .yaml(yaml, c("summary", "metrics"))
    if (is.null(metrics)) {
        return(NULL)
    }
    characters <- metrics[, c("description",
                              "quality_format",
                              "sequence_length")]
    max_size <- metrics[["sequence_length"]] %>%
        gsub(".*-", "", .) %>%
        as.numeric() %>%
        .[] / 10L
    metrics[["library_size"]] <- round(max_size) * 10L
    rownames(metrics) <- metrics[["description"]]
    numerics <- metrics[, setdiff(colnames(metrics), colnames(characters))] %>%
        as.data.frame %>%
        mutate_all(as.numeric)
    bind_cols(characters, numerics) %>%
        as.data.frame %>%
        column_to_rownames("description") %>%
        .[, sort(colnames(.))] %>%
        bind_cols(., characters[, "description", drop = FALSE]) %>%
        DataFrame
}
