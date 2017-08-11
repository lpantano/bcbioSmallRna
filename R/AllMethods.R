#' Count matrix accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `normalized`
#' argument.
#'
#' @rdname counts
#' @docType methods
#'
#' @author Lorena Pantano
#'
#' @param object Object.
#' @param ... Additional parameters.
#'
#' @return Counts matrix
#' @importMethodsFrom DESeq2 counts
#' @export
setMethod("counts", signature("bcbioSmallRnaDataSet"), function(object) {
    counts <- assays(object)[[1]]
    counts
})


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
    cbind(colData(object), metrics) %>% as.data.frame
})


# mirna count data information
#
# @rdname mirna
# @docType methods
#
# @param object [bcbioSmallRnaDataSet] object.
#
# setMethod("mirna", "bcbioSmallRnaDataSet", function(object) {
#     experiments(object)[[1]]
# })

#' bcbioSmallRnaDataSet caller accessors
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
setMethod("bcbio", signature("bcbioSmallRnaDataSet"), function(object, type = "counts") {
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
