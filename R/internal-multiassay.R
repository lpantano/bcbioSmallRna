#' Create multiassay object
#'
#' @rdname multiassay_experiment
#' @noRd
#' @keywords internal
#'
#' @author Lorena Pantano
#'
#'
.multiassay_experiment <- function(exprs, coldata, metadata){
    MultiAssayExperiment(experiments = exprs,
                         colData = coldata,
                         metadata = metadata)
}
