#' Generates WormTensor object
#' A WormTensor object is generated from distance matrices.
#' @param Ds A list containing distance matrices
#' @return An object containing distance matrices and metadata
#' @examples
#' #### test####
#' temp_dl_path <- tempdir()
#' tempfile1 <- file.path(temp_dl_path, "Ds.RData")
#' download.file("https://figshare.com/ndownloader/files/35963780", tempfile1)
#'
#' print(temp_dl_path)
#' dir.exists(temp_dl_path)
#' list.files(temp_dl_path, full.names=TRUE)
#' list.files(temp_dl_path, full.names=TRUE, all.files=TRUE)
#'
#' print(tempfile1)
#' file.exists(tempfile1)
#'
#' tempfile1_back <- paste0(temp_dl_path, "\\Ds.RData")
#' print(tempfile1_back)
#' cat(tempfile1_back)
#' file.exists(tempfile1_back)
#'
#' tempfile1_back2 <- paste0(temp_dl_path, "\\\\Ds.RData")
#' print(tempfile1_back2)
#' cat(tempfile1_back2)
#' file.exists(tempfile1_back2)
#' #### test####
#' worm_download("mSBD", qc = "PASS")$Ds |> as_worm_tensor() -> object
#' @import rTensor
#' @importFrom methods new
#' @importFrom methods is
#' @export
as_worm_tensor <- function(Ds) {
    .check_as_worm_tensor(Ds)
    object <- new("WormTensor")
    object@dist_matrices <- Ds
    object@n_animals <- length(Ds)
    object@union_cellnames <- .union_cellnames(Ds)
    object@n_union_cells <- length(object@union_cellnames)
    object
}

.union_cellnames <- function(Ds) {
    Ds |>
        lapply(function(x) {
            attr(x, "Labels")
        }) |>
        unlist() |>
        unique() |>
        sort() -> res
    res
}

.check_as_worm_tensor <- function(Ds) {
    Ds |>
        lapply(function(x) {
            stopifnot(is(x) == "dist")
        })
}
