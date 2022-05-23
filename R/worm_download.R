#' Downloads distance matrices
#' 28 animals' data including 24 normal and 4 noisy are retrieved from at figshare.
#' @param dist # 追記する
#'
#' @return
#' @export
#'
#' @examples
worm_download <- function(dist="mSBD"){
    Ds <- NULL
    tf <- tempfile()
    if(dist == "mSBD"){
        download.file("https://figshare.com/ndownloader/files/34737937", tf)
    }else if(dist == "Euclid"){
        download.file("https://figshare.com/ndownloader/files/34737931", tf)
    }else{
        stop("Please specify dist as 'mSBD' or 'Euclid'!")
    }
    load(tf)
    Ds
}
