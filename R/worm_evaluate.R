# ここにroxygen2のコメントを書く（下のlibraryは削除）
library("clusterSim") # Pseudo-F
library("clValid") # Connectivity
library("aricode") # ARI
library("cluster") # Silhouette
setMethod("worm_evaluate", "WormTensor",
    function(object, labels){
        # Argument Check
        .check_worm_evaluate(object, labels)
        cluster <- object@clustering
        if(object@clustering_algorithm %in% c("MCMI", "OINDSCAL")){
            data <- object@factor
        }
        if(object@clustering_algorithm == "CSPA"){
            data <- 1 - object@consensus
        }
        #
        # Consistencyをここに追加する
        #
        # Internal Validity Indices（Silhouette係数を追加する）
        psf=.pseudoF(data, cluster)
        cty=.connectivity(data, cluster)
        int_out <- list(PseudoF=psf, Connectivity=cty)
        # External Validity Indices
        if(!is.null(labels)){
            ext_out <- lapply(labels, function(l){
                list(
                    Fmeasure=.fmeasure(cluster, l),
                    Entropy=.entropy(cluster, l),
                    Purity=.purity(cluster, l),
                    ARI=ARI(cluster, l))
            })
            names(ext_out) <- names(labels)
        }else{
            ext_out <- list(Fmeasure=NULL, Entropy=NULL, Purity=NULL)
        }
        # Ouput
        out <- list(internal=int_out, external=ext_out)
        object@eval <- out
        object
    }
)

.check_worm_evaluate <- function(object, labels){
    # Backword Check
    if(prod(object@clustering) == 1){
        msg <- paste0("Perform worm_clustering() first.")
        message(msg)
    }
    # Argument Check
    if(!is.null(labels)){
        if(!is.list(labels)){
            stop("Specify labels as a list.")
        }
    }
}

######### Internal Validity Indices (w/o Labels) #########
.pseudoF <- function(data, cluster){
    index.G1(data, cluster)
}

.connectivity <- function(data, cluster){
    # Distance matrix
    Dist <- dist(data, method="euclidean")
    # Connectivity
    connectivity(Dist, cluster)
}

######### External Validity Indices (w Labels) #########
.fmeasure <- function(cluster, label){
    ctbl <- table(cluster, label)
    # All combination of Recall
    R <- ctbl / colSums(ctbl)
    # All combination of Precision
    P <- ctbl / rowSums(ctbl)
    # All combination of F-measure
    F <- 2 * R * P / (R + P)
    # NaN => 0
    F[is.nan(F)] <- 0
    # Weight
    w <- apply(ctbl, 2, sum) / sum(ctbl)
    # Total Micro-averaged F value
    sum(w * apply(F, 2, max))
}

.entropy <- function(cluster, label){
    # Cross tabulation
    ctbl <- table(cluster, label)
    # Weight
    w <- apply(ctbl, 1, sum) / sum(ctbl)
    # Total Entropy
    sum(w * apply(ctbl, 1, .calcEntropy0))
}

.calcEntropy0 <- function(pv){
    p1 <- pv / sum(pv)
    p2 <- p1[p1 !=0]
    # Entropy each column of ctbl
    - sum(p2 * log2(p2))
}

.purity <- function(cluster, label){
    # Cross tabulation
    ctbl <- table(cluster, label)
    # Weight
    w <- apply(ctbl, 1, sum) / sum(ctbl)
    # Purity
    sum(w * apply(ctbl, 1, max) / rowSums(ctbl))
}
