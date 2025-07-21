#' Downloads distance matrices
#' 28 animals' data including 24 normal and 4 noisy are retrieved from figshare.
#' @param distance "mSBD" or "Euclid" can be specified
#' @param qc "PASS" or "WARN" or "FAIL" can be specified. "PASS" downloads
#' 24 data except 4 noisy data. "WARN" downloads 27 data except 1 noisy data.
#' "FAIL" downloads all 28 data.
#' @return A List of containing distance matrices. The list also includes
#' metadata for each animals.
#' @examples
#' \donttest{
#'     Ds_Euclid <- worm_download("Euclid", qc = "WARN")
#'     Ds_mSBD <- worm_download("mSBD", qc = "PASS")
#' }
#' @importFrom usedist dist_subset
#' @importFrom utils download.file
#' @importFrom utils read.csv
#' @importFrom curl has_internet
#' @export
worm_download <- function(distance = c("mSBD", "Euclid"),
                          qc       = c("PASS", "WARN", "FAIL")) {
  distance <- match.arg(distance)
  qc       <- match.arg(qc)
  
  # CRAN fails internet access -> stop gracefully
  if (!curl::has_internet()) {
    stop(
      "WormTensor requires internet to fetch data, but no connection was detected.\n",
      "Please check your network or download data manually from:\n",
      "  https://figshare.com/articles/…/35963780 or 35963777"
    )
  }
  
  safe_download <- function(url, dest, mode = NULL) {
    tryCatch({
      if (is.null(mode)) {
        utils::download.file(url, dest)
      } else {
        utils::download.file(url, dest, mode = mode)
      }
    }, error = function(e) {
      stop(
        "Failed to download WormTensor data from:\n  ", url, "\n",
        "Please check your internet connection or visit the URL manually."
      )
    })
  }
  
  tmpdir    <- tempdir()
  file_rda  <- file.path(tmpdir, "Ds.RData")
  file_lbl  <- file.path(tmpdir, "labels.csv")
  
  # Download
  url_data <- switch(distance,
                     mSBD   = "https://figshare.com/ndownloader/files/35963780",
                     Euclid = "https://figshare.com/ndownloader/files/35963777"
  )
  safe_download(
    url_data, file_rda,
    mode = if (.Platform$OS.type == "windows") "wb" else NULL
  )
  load(file_rda)
  
  # QC filtering
  idx <- switch(qc,
                PASS = which(.qcresult %in% "PASS"),
                WARN = which(.qcresult %in% c("PASS","WARN")),
                FAIL = seq_along(.qcresult)
  )
  
  # Labels
  safe_download(
    "https://figshare.com/ndownloader/files/36186483",
    file_lbl,
    mode = if (.Platform$OS.type == "windows") "wb" else NULL
  )
  labels <- utils::read.csv(file_lbl)
  
  ## Output
  Ds_f <- lapply(Ds, .filter_cellnames)
  list(Ds     = Ds_f[idx],
       labels = labels)
}

.qcresult <- rep("", 28)
.qcresult[c(1:2, 4:7, 9:19, 21:24, 26:28)] <- "PASS"
.qcresult[c(3, 8, 25)] <- "WARN"
.qcresult[c(20)] <- "FAIL"
.filter_cellnames <- function(X) {
    D <- X
    D_cell <- attr(D, "Labels")
    D_cell_f <- D_cell[grep("^[0-9]", D_cell, invert = TRUE)]
    D_f <- dist_subset(D, D_cell_f)
    D_f
}
