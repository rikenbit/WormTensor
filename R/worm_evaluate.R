#' Title
#'
#' @param WormTensor
#'
#' @return
#' @examples
#' # Pipe Operation
#' worm_download("mSBD", qc="PASS")$Ds |>
#'     as_worm_tensor() |>
#'         worm_membership(k=6) |>
#'             worm_clustering() -> object
#' worm_evaluate(object) -> object_internal
#' labels <- list(
#'     label1 = sample(3, length(object@clustering), replace=TRUE),
#'     label2 = sample(4, length(object@clustering), replace=TRUE),
#'     label3 = sample(5, length(object@clustering), replace=TRUE))
#' worm_evaluate(object, labels) -> object_external
#' Ds_mSBD <- worm_download("mSBD", qc="PASS")
#' # label1 is used Consistency
#' labels <- list(
#'     label1 = replace(Ds_mSBD$labels$Class, which(is.na(Ds_mSBD$labels$Class)), "NA"),
#'     label2 = sample(4, length(object@clustering), replace=TRUE),
#'     label3 = sample(5, length(object@clustering), replace=TRUE))
#' worm_evaluate(object, labels) -> object_external_Class
#' @importFrom clusterSim index.G1 # Pseudo-F
#' @importFrom clValid connectivity # Connectivity
#' @importFrom aricode ARI # ARI
#' @importFrom cluster silhouette # Silhouette
#' @export
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
        # Consistencyの計算で使う個体ごとのクラスタリング結果
        Cs <- lapply(object@dist_matrices, function(d, k){
            cutree(hclust(d, method="ward.D2"), k)
        }, k=object@k)
        # labelsの1つめのリストのConsistencyを取得
        consistency=.consistency(object, labels, Cs)[[1]]$Consistency

        # No. of identified cells（各個体のCellCountを追加）
        no_identified=.no_identified(object)

        # Internal Validity Indices（Silhouette係数を追加する）
        silhouette=.silhouette_cell(object)

        cellwise <- list(consistency=consistency,
                         no_identified=no_identified,
                         silhouette=silhouette)

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
        out <- list(internal=int_out, external=ext_out, cellwise=cellwise)
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
######### consistency #########
.consistency <- function(object, labels, Cs){
    consistency_l <- lapply(labels, function(l){
        df_eval_label <- data.frame(CellType = object@union_cellnames,
                                    Classes = l,
                                    stringsAsFactors = FALSE)

        df_count_union_list <- lapply(Cs, function(cs){
            df_cls <- data.frame(CellType = names(cs),
                                 Cluster = cs,
                                 stringsAsFactors = FALSE,
                                 row.names = NULL)
            df_cls_label <- merge(df_cls,
                                  df_eval_label,
                                  by.x = "CellType",
                                  by.y = "CellType",
                                  all.x = TRUE)
            label <- df_cls_label$Classes
            cluster <- df_cls_label$Cluster
            H_label <- .H(label)
            H_cluster <- .H(cluster)

            other_member <- .other_member(H_label)

            # カウント
            count1 <- rep(0, length=length(label))
            for(i in seq_len(nrow(H_label))){
                idx1 <- which(H_label[i, ] == 1)
                idx2 <- which(H_cluster[i, ] == 1)
                count1[i] <- sum(H_label[,idx1] * H_cluster[,idx2]) - 1
            }
            # 正規化
            count1 <- count1 / other_member
            # 0/0のNaN対策
            count1 <- ifelse(is.nan(count1), 0, count1)
            # 1個体分 count1
            df_count <- cbind(df_cls_label, Count=count1)
            # object@union_cellnamesにcount1を埋める
            df_union <- data.frame(CellType = object@union_cellnames,
                                   stringsAsFactors = FALSE)
            df_count_union <- merge(df_union,
                                    df_count,
                                    by.x = "CellType",
                                    by.y = "CellType",
                                    all.x = TRUE)
            df_count_union_ <- df_count_union[,c("Count")]
            df_count_union_[is.na(df_count_union_)] <- 0
            return(df_count_union_)
        })
        consistency <- Reduce("+",df_count_union_list)
        df_count_sum <- data.frame(CellType = object@union_cellnames,
                                   Consistency = consistency,
                                   stringsAsFactors = FALSE)
        return(df_count_sum)
    })
    return(consistency_l)
}
# 一回、接続行列に変換
.H <- function(vec){
    uniq_vec <- unique(vec)
    out <- matrix(0, nrow=length(vec), ncol=length(uniq_vec))
    rownames(out) <- names(vec)
    colnames(out) <- uniq_vec
    for(i in seq_len(length(vec))){
        idx <- which(vec[i] == uniq_vec)
        out[i, idx] <- 1
    }
    out
}
# 自分以外のメンバー数 正規化用
.other_member <- function(H_label){
    out <- rep(0, length=nrow(H_label))
    for(i in seq_len(nrow(H_label))){
        idx1 <- which(H_label[i, ] == 1)
        out[i] <- sum(H_label[,idx1]) - 1
    }
    out
}
######### no_identified #########
.no_identified <- function(object){
    lapply(object@dist_matrices, function(d){
        attr(d, "Labels")
    }) |>
        unlist() -> all_cellname
    all_cellname |>
        table() |>
            as.numeric() -> cell_count
    return(cell_count)
}
######### no_identified #########
.silhouette_cell <- function(object){
    algorithm <- object@clustering_algorithm
    sil_dist <- .distFunc[[algorithm]]
    cls <- object@clustering

    sil <- silhouette(cls, sil_dist)
    rownames(sil) <- object@union_cellnames
    return(sil)
}
# Selecting how to create a distance matrix object
.distFunc <- list(
    "CSPA" = as.dist(1 - object@consensus),
    "OINDSCAL" = dist(object@factor),
    "MCMI" = dist(object@factor)
)
