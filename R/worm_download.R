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
