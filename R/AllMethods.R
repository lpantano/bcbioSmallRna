#' Sample metrics accesor
#'
#' @rdname metrics
#' @docType methods
#'
#' @param object [bcbioSmallRnaDataSet] object.
#' @export
setMethod("metrics", signature("bcbioSmallRnaDataSet"), function(object) {
    metrics <- metadata(object)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    right_join(as.data.frame(colData(object)),
               as.data.frame(metrics), by = c("sample" = "description")) %>% as.data.frame
})


#' mirna count data information
#
#' @rdname mirna
#' @docType methods
#' @aliases cluster trna isomir
#' @param object [bcbioSmallRnaDataSet] object.
#' @param type Type of count data to retrieve. (raw rlog coldata)
#' @export
setMethod("mirna", "bcbioSmallRnaDataSet", function(object, type="raw") {
    if (type == "coldata")
        return(colData(experiments(object)[["mirna"]]))
    if (type %in% names(assays(experiments(object)[["mirna"]])))
        assays(experiments(object)[["mirna"]])[[type]]
    else stop(paste(type, "not found."))
})
#' @rdname mirna
#' @export
setMethod("isomir", "bcbioSmallRnaDataSet", function(object, type="raw") {
    if (type == "coldata")
        return(colData(experiments(object)[["isomirs"]]))
    if (type %in% names(assays(experiments(object)[["isomirs"]])))
        assays(experiments(object)[["isomirs"]])[[type]]
    else stop(paste(type, "not found."))
})
#' @rdname mirna
#' @export
setMethod("cluster", "bcbioSmallRnaDataSet", function(object, type="raw") {
    if (type %in% names(assays(experiments(object)[["cluster"]])))
        assays(experiments(object)[["cluster"]])[[type]]
    else stop(paste(type, "not found"))
})
#' @rdname mirna
#' @export
setMethod("trna", "bcbioSmallRnaDataSet", function(object, type="raw") {
    if (type %in% names(assays(experiments(object)[["trna"]])))
        assays(experiments(object)[["trna"]])[[type]]
    else stop(paste(type, "not found"))
})

#' `bcbioSmallRnaDataSet` caller accessors
#'
#'
#' Additional obejcts of interest include:
#'
#' - [isomiRs](http://bioconductor.org/packages/3.6/bioc/html/isomiRs.html)
#'
#' @rdname bcbio
#' @docType methods
#'
#' @param object Object.
#' @param type Type of count data to retrieve.
#' @param value An integer matrix or other object.
#' @param ... Additional arguments.
#'
#' @return caller object.
#' @export
setMethod("bcbio", signature("bcbioSmallRnaDataSet"),
          function(object, type = "counts") {
    if (type %in% names(slot(object, "callers"))) {
        slot(object, "callers")[[type]]
    } else {
        stop(paste(type, "not found"))
    }
})
#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(object = "bcbioSmallRnaDataSet", value = "ANY"),
    function(object, type = "counts", value) {

        slot(object, "callers")[[type]] <- value
        object
    })



#' Sample adapter accesor
#'
#' @rdname adapter
#' @docType methods
#'
#' @param object [bcbioSmallRnaDataSet] object.
#' @param value A [data.frame] with adapter removal informaion.
#' @param ... Additional arguments.
#'
#' @return [list] or [bcbioSmallRnaDataSet].
#' @export
setMethod("adapter", signature("bcbioSmallRnaDataSet"), function(object) {
    adapter <- metadata(object)[["adapter"]]
    if (is.null(adapter)) return(NULL)
    adapter
})
#' @rdname adapter
#' @export
setMethod(
    "adapter<-",
    signature(object = "bcbioSmallRnaDataSet", value = "list"),
    function(object, value) {
        metadata(object)[["adapter"]] <- value
        object
    })
