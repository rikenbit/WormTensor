#' Generates WormTensor object
#' A WormTensor object is generated from distance matrices.
#' @param Ds A List of containing distance matrices. The list also includes
#' metadata for each animals.
#' @return An object containing distance matrices and metadata
#' @examples
#' worm_download("mSBD", qc="PASS")$Ds |> as_worm_tensor() -> object
#' @import rTensor
#' @importFrom methods new
#' @export
as_worm_tensor <- function(Ds){
    .check_as_worm_tensor(Ds)
    object <- new("WormTensor")
    object@dist_matrices <- Ds
    object@n_animals <- length(Ds)
    object@union_cellnames <- .union_cellnames(Ds)
    object@n_union_cells <- length(object@union_cellnames)
    object
}

.union_cellnames <- function(Ds){
    # for worm_download data
    sort(unique(unlist(lapply(Ds, function(x){attr(x, "Labels")})))) -> res
    # for sample data (worm_distance.R's Toy data)
    if(is.null(res)){
        sort(unique(unlist(lapply(Ds, function(x){attr(x, "dimnames")})))) -> res
    }
    res
}

.check_as_worm_tensor <- function(Ds){
    lapply(Ds, function(x){
        stopifnot(is(x) == "dist")
    })
}
