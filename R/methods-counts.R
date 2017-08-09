#' Count matrix accessors
#'
#' By default, [counts()] returns the raw counts. Normalized counts, including
#' transcripts per million (TPM) can be accessed using the `normalized`
#' argument.
#'
#' @rdname counts
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param ... Additional parameters.
#'
#' @return Counts matrix
#' @export
#'
#' @examples
#' data(dummy)
#' # Raw counts
#' ma <- counts(bcb)
setMethod("counts", "bcbioSmallRnaDataSet", function(object, normalized = FALSE) {
    counts <- assays(object)[[1]]

    counts
})
