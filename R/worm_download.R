#' Downloads distance matrices
#' 28 animals' data including 24 normal and 4 noisy are retrieved from at figshare.
#' @param dist # 追記する
#'
#' @return
#' @export
#'
#' @examples
worm_download <- function(distance=c("mSBD", "Euclid"),
                          qc=c("PASS", "WARN", "FAIL")){
    # Argument Check
    distance <- match.arg(distance)
    qc <- match.arg(qc)
    # Distance matrices
    Ds <- NULL
    tmpdir <- tempdir()
    tempfile1 <- paste0(tmpdir, "/Ds.RData")
    if(distance == "mSBD"){
        download.file(
            "https://figshare.com/ndownloader/files/34737937",
            tempfile1)
    }else if(distance == "Euclid"){
        download.file(
            "https://figshare.com/ndownloader/files/34737931",
            tempfile1)
    }else{
        stop("Please specify distance as 'mSBD' or 'Euclid'!")
    }
    load(tempfile1)
    if(qc == "PASS"){
        idx <- which(.qcresult %in% c("PASS"))
    }
    if(qc == "WARN"){
        idx <- which(.qcresult %in% c("PASS", "WARN"))
    }
    if(qc == "FAIL"){
        idx <- seq(28)
    }
    # Labels
    tempfile2 <- paste0(tmpdir, "/labels.csv")
    download.file("https://figshare.com/ndownloader/files/34398068",
                  tempfile2)
    labels <- read.csv(tempfile2)
    # Output
    list(Ds=Ds[idx], labels=labels)
}

.qcresult <- rep("", 28)
.qcresult[c(1:2,4:7,9:19,21:24,26:28)] <- "PASS"
.qcresult[c(3,8,25)] <- "WARN"
.qcresult[c(20)] <- "FAIL"
