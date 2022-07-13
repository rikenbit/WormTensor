#' Title
#'
#' @param WormTensor
#'
#' @return
#' @examples
#' # Temporary directory to save figures
#' out.dir <- tempdir()
#'
#' #evaluate
#' worm_download("mSBD", qc="PASS")$Ds |>
#'     as_worm_tensor() |>
#'         worm_membership(k=6) |>
#'             worm_clustering() -> object
#' Ds_mSBD <- worm_download("mSBD", qc="PASS")
#' labels <- list(
#'     label1 = replace(Ds_mSBD$labels$Class, which(is.na(Ds_mSBD$labels$Class)), "NA"),
#'     label2 = sample(4, length(object@clustering), replace=TRUE),
#'     label3 = sample(5, length(object@clustering), replace=TRUE))
#'
#' # Pipe Operation
# worm_download("mSBD", qc="PASS")$Ds |>
#     as_worm_tensor() |>
#     worm_membership(k=6) |>
#     worm_clustering() |>
#     worm_evaluate(labels) |>
#     worm_visualize("tSNE",out.dir) -> object
#' @importFrom
#' @import ggplot2
#' @import Rtsne
#' @import uwot
#' @import factoextra
#' @export
setMethod("worm_visualize", "WormTensor",
    function(object, algorithm, out.dir){
    # Argument Check
    algorithm <- match.arg(algorithm)
    .check_worm_visualize(object, out.dir)
    # data forDimensional Reduction
    if(object@clustering_algorithm %in% c("MCMI", "OINDSCAL")){
        data <- object@factor
    }
    if(object@clustering_algorithm == "CSPA"){
        data <- 1 - object@consensus
    }
    # dist for Dimensional Reduction
    if(object@clustering_algorithm %in% c("MCMI", "OINDSCAL")){
        cls_dist <- dist(data)
    }
    if(object@clustering_algorithm == "CSPA"){
        cls_dist <- as.dist(data)
    }
    # Random number fixing
    set.seed(1234)
    if(algorithm == "tSNE"){
        twoD <- Rtsne(cls_dist,
                      is_distance=TRUE,
                      # perplexity = 15,
                      # max_iter = 1000,
                      check_duplicates=FALSE
                      )$Y
    }
    if(algorithm == "UMAP"){
        # cf. https://github.com/jlmelville/uwot/issues/22
        twoD <- uwot::umap(X = NULL,
                           metric = "precomputed",
                           # n_neighbors = 15,
                           # n_components = 2,
                           nn_method = uwot:::dist_nn(cls_dist,
                                                      k = attr(cls_dist, "Size")),
                           )
    }

    # ここに画像をout.dir以下に出力するコードを書いていく
    # Make figures dir
    if(!dir.exists(paste0(out.dir, "/figures"))){
        dir.create(paste0(out.dir, "/figures"))
    }
    #### 1. 細胞ごとのシルエット図（例: 論文 Figure 2）####
    sil <- object_eval@eval$cellwise$silhouette
    gg_sil <- fviz_silhouette(sil) +
        labs(y = "Silhouette width",
             x = "",
             title = "",
             color ="Cluster",
             fill = "Cluster") +
        theme(text = element_text(size = 90))
    # save (silhouette plot of each cell)
    ggsave(filename = paste0(out.dir, "/figures/silhouette.png"),
           plot = gg_sil,
           dpi = 100,
           width = 30.0,
           height = 20.0,
           limitsize = FALSE)

    ##### 2. 次元圧縮図に色を反映させたもの（例: 論文 Figure 3,4）####


    # 3. 重み/ARIと同定細胞数の関係（例: 論文 Figure 6a）

    # Output
    object
    }
)

.check_worm_visualize <- function(object, out.dir){
    # Backword Check
    if(length(object@eval) == 0){
        stop("Perform worm_eval first.")
    }
    # Argument Check
    if(!is.null(out.dir)){
        if(!file.exists(out.dir)){
            stop("Specify a valid out.dir.")
        }
    }
}
