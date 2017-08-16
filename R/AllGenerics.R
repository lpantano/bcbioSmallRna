#' @rdname adapter
#' @usage NULL
#' @export
setGeneric("adapter", function(object, ...) standardGeneric("adapter"))
#' @rdname adapter
#' @usage NULL
#' @export
setGeneric("adapter<-", function(object, ..., value)
    standardGeneric("adapter<-"))



#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))
#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



#' @rdname metrics
#' @usage NULL
#' @export
setGeneric("metrics", function(object) standardGeneric("metrics"))



#' @rdname mirna
#' @usage NULL
#' @export
setGeneric("mirna", function(object, ...) standardGeneric("mirna"))
#' @rdname mirna
#' @usage NULL
#' @export
setGeneric("cluster", function(object, ...) standardGeneric("cluster"))
#' @rdname mirna
#' @usage NULL
#' @export
setGeneric("trna", function(object, ...) standardGeneric("trna"))
#' @rdname mirna
#' @usage NULL
#' @export
setGeneric("isomir", function(object, ...) standardGeneric("isomir"))