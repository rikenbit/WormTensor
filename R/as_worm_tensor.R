#' Generates WormTensor object
#' A WormTensor object is generated from distance matrices.
#' @param Ds A list containing distance matrices
#' @return An object containing distance matrices and metadata
#' @examples
#' #### test####
#' temp_dl_path <- tempdir()
#' print(temp_dl_path)
#' dir.exists(temp_dl_path)
#'
#' tempfile1 <- file.path(temp_dl_path, "Ds.RData")
#' download.file("https://figshare.com/ndownloader/files/35963780", tempfile1, mode="wb")
#'
#' print(tempfile1)
#' file.exists(tempfile1)
#' file.access(tempfile1,mode=0)
#' file.access(tempfile1,mode=1)
#' file.access(tempfile1,mode=2)
#' file.access(tempfile1,mode=4)
#'
#' options()$download.file.method
#' file.info(tempfile1)
#'
#' load(tempfile1)
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
