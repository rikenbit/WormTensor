# ここにroxygen2のコメントを書く（下のlibraryは削除）
library("ggplot2")
library("Rtsne")
library("uwot")
setMethod("worm_visualize", "WormTensor",
    function(object, algorithm, out.dir){
    # Argument Check
    algorithm <- match.arg(algorithm)
    .check_worm_visualize(object, out.dir)
    # Dimensional Reduction
    if(object@clustering_algorithm %in% c("MCMI", "OINDSCAL")){
        data <- object@factor
    }
    if(object@clustering_algorithm == "CSPA"){
        data <- as.dist(1 - object@consensus)
    }
    if(algorithm == "tSNE"){
        if("dist" %in% is(data)){
            twoD <- Rtsne(data, is_distance=TRUE, check_duplicates=FALSE)$Y
        }else{
            twoD <- Rtsne(data, check_duplicates=FALSE)$Y
        }
    }
    if(algorithm == "UMAP"){
        twoD <- umap(data)
    }
    # ここに画像をout.dir以下に出力するコードを書いていく
    # 1. 細胞ごとのシルエット図（例: 論文 Figure 2）
    #   （object_internal@eval$internal$Silhouetteの値を利用）

    # 2. 次元圧縮図に色を反映させたもの（例: 論文 Figure 3,4）

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
