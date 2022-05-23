# ここにroxygen2のコメントを書く（下のlibraryは削除）
library("rTensor")
as_worm_tensor <- function(Ds){
    object <- new("WormTensor")
    object@dist_matrices <- Ds
    object@n_animals <- length(Ds)
    object@union_cellnames <- .union_cellnames(Ds)
    object@n_union_cells <- length(object@union_cellnames)
    object
}

.union_cellnames <- function(Ds){
    sort(unique(unlist(lapply(Ds, names))))
}
