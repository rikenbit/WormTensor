#' Generates WormTensor object
#' A WormTensor object is generated from distance matrices.
#' @param Ds A List of containing distance matrices. The list also includes
#' metadata for each animals.
#' @return An object containing distance matrices and metadata
#' @examples
#' worm_download("Euclid", qc="WARN")$Ds |> as_worm_tensor() -> object
#' @import rTensor
#' @export
as_worm_tensor <- function(Ds){
    object <- new("WormTensor")
    object@dist_matrices <- Ds
    object@n_animals <- length(Ds)
    object@union_cellnames <- .union_cellnames(Ds)
    object@n_union_cells <- length(object@union_cellnames)
    object
}

.union_cellnames <- function(Ds){
    sort(unique(unlist(lapply(Ds, function(x){attr(x, "Labels")}))))
}
