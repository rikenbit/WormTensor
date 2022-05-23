# ここにroxygen2のコメントを書く（下のlibraryは削除）
library("rTensor")
setClass("WormTensor",
    slots = c(
        # Filled by as_worm_tensor()
        dist_matrices = "list",
        n_animals = "numeric",
        union_cellnames = "character",
        n_union_cells = "numeric",
        # Filled by worm_membership()
        membership_tensor = "Tensor",
        k = "numeric",
        # Filled by worm_clustering()
        clustering_algorithm = "character",
        clustering = "numeric",
        weight = "numeric",
        factor = "matrix",
        consensus = "matrix",
        # Filled by worm_evaluate()
        eval = "list",
        # Filled by worm_visualize()
        dimension_reduction_algorithm = "character"))
