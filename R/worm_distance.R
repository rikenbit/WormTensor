# ここにroxygen2のコメントを書く（下のlibraryは削除）
library("dtwclust")
worm_distance <- function(data, distance=c("mSBD", "SBD", "Euclid")){
    # Argument Check
    .check_worm_distance(data)
    distance <- match.arg(distance)
    # Calculate Distance
    Ds <- lapply(data, .distFunc[[distance]])
    # Output
    Ds
}

.check_worm_distance <- function(data){
    stopifnot(is.list(data))
    if(length(data) <= 2){
        stop("The number of matrices are too small!")
    }
    stopifnot(all(unlist(lapply(data, is.numeric))))
}

.SBDMatrix <- function(X){
    indices  <- t(combn(seq(nrow(X)), 2))
    indices <- indices[,2:1]
    distances <- apply(indices, 1, function(xx){
        SBD(X[xx[1],], X[xx[2],])$dist
    })
    out <- matrix(0, nrow=nrow(X), ncol=nrow(X))
    out[indices] <- distances
    as.dist(out)
}

.mSBDMatrix <- function(X){
    indices  <- t(combn(seq(nrow(X)), 2))
    indices <- indices[,2:1]
    distances <- apply(indices, 1, function(xx){
        .mSBD(X[xx[1],], X[xx[2],])$dist
    })
    out <- matrix(0, nrow=nrow(X), ncol=nrow(X))
    out[indices] <- distances
    as.dist(out)
}

.mSBD <- function(x, y, znorm = FALSE, error.check = TRUE,
    return.shifted = TRUE){
    nx <- length(x)
    ny <- length(y)
    if(nx > ny){
        swap <- x
        x <- y
        y <- swap
    }else{
        swap <- NULL
    }
    if(znorm){
        CCseq <- NCCc(zscore(x, error.check = FALSE),
                      zscore(y, error.check = FALSE),
                      error.check = FALSE)
    }else{
        CCseq <- NCCc(x, y, error.check = FALSE)
    }
    m <- max(abs(CCseq))
    if(!return.shifted){
        return(1 - m) # nocov
    }
    shift <- which.max(abs(CCseq)) - max(nx, ny)
    if(is.null(swap)){
        if(shift < 0L){
            yshift <- y[(-shift + 1L):ny]
        }else{
            yshift <- c(rep(0, shift), y)
        }
    }else{
        if(shift < 0L){
            yshift <- c( rep(0, -shift), x )
        }else{
            yshift <- x[(shift + 1L):ny]
        }
    }
    nys <- length(yshift)
    if (nys < nx){
        yshift <- c( yshift, rep(0, nx-nys) )
    }else{
        yshift <- yshift[1L:nx]
    }
    # return
    list(dist = 1 - m, yshift = yshift)
}

.distFunc <- list(
    "Euclid" = dist,
    "SBD" = .SBDMatrix,
    "mSBD" = .mSBDMatrix
)
