#' Sample metrics accesor
#'
#' @rdname metrics
#' @docType methods
#'
#' @param object [bcbioSmallRnaDataSet] object.
#'
#' @export
setMethod("metrics", "bcbioSmallRnaDataSet", function(object) {
    metrics <- metadata(object)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    cbind(colData(object), metrics) %>% as.data.frame
})
