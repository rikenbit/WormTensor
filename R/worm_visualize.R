#' Plots evaluation result
#' A visualization result is generated from a WormTensor object
#' @param object WormTensor object with a result of worm_evaluate
#' @param algorithm Dimensional reduction methods
#' @param out.dir Output directory (default: tempdir())
#' @param seed Arguments passed to set.seed (default: 1234)
#' @param tsne.dims Output dimensionality (default: 2)
#' @param tsne.perplexity Perplexity paramete (default: 15)
#' @param tsne.verbose logical; Whether progress updates should be printed
#' (default: TRUE)
#' @param tsne.max_iter Number of iterations (default: 1000)
#' @param umap.n_neighbors The size of local neighborhood (default: 15)
#' @param umap.n_components The dimension of the space to embed into
#' (default: 2)
#' @param silhouette.summary logical; If true a summary of
#' cluster silhouettes are printed.
#' @return Silhouette plots. ARI with a merge result and each animal(with MCMI).
#'  Dimensional reduction Plots colored by cluster, no. of identified cells,
#'  consistency(with labels), Class_label(with labels).
#' @examples
#' # Temporary directory to save figures
#' out.dir <- tempdir()
#'
#' # labels
#' worm_download("mSBD", qc="PASS")$Ds |>
#'     as_worm_tensor() |>
#'         worm_membership(k=6) |>
#'             worm_clustering() -> object
#' Ds_mSBD <- worm_download("mSBD", qc="PASS")
#' labels <- list(
#'     label1 = replace(Ds_mSBD$labels$Class,
#'     which(is.na(Ds_mSBD$labels$Class)),
#'     "NA"),
#'     label2 = sample(4, length(object@clustering), replace=TRUE),
#'     label3 = sample(5, length(object@clustering), replace=TRUE))
#'
#' # Pipe Operation (no labels)
#' worm_download("mSBD", qc="PASS")$Ds |>
#'    as_worm_tensor() |>
#'    worm_membership(k=6) |>
#'    worm_clustering() |>
#'    worm_evaluate() |>
#'    worm_visualize("tSNE",out.dir) -> object_no_labels
#'
#' # Pipe Operation (with labels)
#' worm_download("mSBD", qc="PASS")$Ds |>
#'    as_worm_tensor() |>
#'    worm_membership(k=6) |>
#'    worm_clustering() |>
#'    worm_evaluate(labels) |>
#'    worm_visualize("tSNE",out.dir) -> object_labels
#' @import ggplot2
#' @importFrom Rtsne Rtsne
#' @import uwot
#' @importFrom factoextra fviz_silhouette
#' @importFrom ggrepel geom_label_repel
#' @importFrom cowplot align_plots
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @importFrom cowplot theme_half_open
#' @export
setMethod("worm_visualize", "WormTensor",
    function(object,
             algorithm,
             out.dir,
             seed,
             tsne.dims,
             tsne.perplexity,
             tsne.verbose,
             tsne.max_iter,
             umap.n_neighbors,
             umap.n_components,
             silhouette.summary){
    # Argument Check
    algorithm <- match.arg(algorithm)
    object@dimension_reduction_algorithm <- algorithm
    .check_worm_visualize(object,
                          out.dir,
                          seed,
                          tsne.dims,
                          tsne.perplexity,
                          tsne.verbose,
                          tsne.max_iter,
                          umap.n_neighbors,
                          umap.n_components,
                          silhouette.summary)
    # data for Dimensional Reduction
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
    set.seed(seed)
    # Dimensional Reduction
    if(algorithm == "tSNE"){
        twoD <- Rtsne(cls_dist,
                      is_distance=TRUE,
                      dims = tsne.dims,
                      perplexity = tsne.perplexity,
                      verbose = tsne.verbose,
                      max_iter = tsne.max_iter,
                      check_duplicates=FALSE
                      )$Y
        cord_x <- c("t-SNE-1")
        cord_y <- c("t-SNE-2")
    }
    if(algorithm == "UMAP"){
        # cf. https://github.com/jlmelville/uwot/issues/22
        twoD <- uwot::umap(X = NULL,
                           metric = "precomputed",
                           n_neighbors = umap.n_neighbors,
                           n_components = umap.n_components,
                           nn_method = uwot:::dist_nn(cls_dist,
                                                      k = attr(cls_dist, "Size")
                                                      )
                           )
        cord_x <- c("UMAP-1")
        cord_y <- c("UMAP-2")
    }

    # Make figures directory
    if(!dir.exists(paste0(out.dir, "/figures"))){
        dir.create(paste0(out.dir, "/figures"))
    }
    #### 1. 細胞ごとのシルエット図（例: 論文 Figure 2）####
    sil <- object@eval$cellwise$silhouette
    gg_sil <- fviz_silhouette(sil, print.summary = silhouette.summary) +
        labs(y = "Silhouette width",
             x = "",
             # title = "",
             color ="Cluster",
             fill = "Cluster") +
        theme(text = element_text(size = 90))
    # save (silhouette plot of each cell)
    ggsave(filename = paste0(out.dir, "/figures/Silhouette.png"),
           plot = gg_sil,
           dpi = 100,
           width = 30.0,
           height = 20.0,
           limitsize = FALSE)

    ##### 2. 次元圧縮図に色を反映させたもの（例: 論文 Figure 3,4）####
    ##### Cluster#####
    df_cls <- data.frame(cord_1 = twoD[,1],
                         cord_2 = twoD[,2],
                         cell_type = object@union_cellnames,
                         Cluster = object@clustering,
                         stringsAsFactors = FALSE
                         )
    gg_cls <- ggplot(df_cls,
                     aes(x = cord_1,
                         y = cord_2,
                         label = cell_type,
                         color = factor(Cluster)
                         )
                     ) +
        labs(color = "Cluster") +
        geom_point(size = 6.0,
                   alpha = 0.6) +
        geom_label_repel(max.overlaps = Inf,
                         min.segment.length = 0,
                         size = 9.0,
                         force = 6.0) +
        theme(text = element_text(size = 60)) +
        labs(x = cord_x,
             y = cord_y)
    ggsave(filename = paste0(out.dir, "/figures/Cluster.png"),
           plot = gg_cls,
           dpi = 100,
           width = 25.0,
           height = 20.0,
           limitsize = FALSE)

    ##### No. of cells#####
    df_cc <- data.frame(cord_1 = twoD[,1],
                        cord_2 = twoD[,2],
                        cell_type = object@union_cellnames,
                        cell_count = object@eval$cellwise$no_identified,
                        stringsAsFactors = FALSE
                        )
    gg_cc <- ggplot(df_cc,
                    aes(x = cord_1,
                        y = cord_2,
                        label = cell_type,
                        color = cell_count
                        )
        ) +
        scale_color_viridis_c(option = "D") +
        labs(color = "No. of cells") +
        geom_point(size = 6.0,
                   alpha = 0.6) +
        geom_label_repel(max.overlaps = Inf,
                         min.segment.length = 0,
                         size = 9.0,
                         force = 6.0) +
        theme(text = element_text(size = 60)) +
        labs(x = cord_x,
             y = cord_y) +
        theme(legend.key.height = unit(1.5, "cm")) +
        theme(legend.key.width = unit(1.5, "cm"))
    ggsave(filename = paste0(out.dir, "/figures/no_identified.png"),
           plot = gg_cc,
           dpi = 100,
           width = 25.0,
           height = 20.0,
           limitsize = FALSE)
    ##### consistency  (plot each label)#####
    if(!is.null(object@eval$cellwise$consistency)){
        con_list <- object@eval$cellwise$consistency
        for (x in seq_along(con_list)) {
            con_name <- names(con_list[x])
            con_value <- con_list[[x]]
            df_con <- data.frame(cord_1 = twoD[,1],
                            cord_2 = twoD[,2],
                            cell_type = object@union_cellnames,
                            consistency = con_value,
                            stringsAsFactors = FALSE
                            )
            gg_con <- ggplot(df_con,
                             aes(x = cord_1,
                                 y = cord_2,
                                 label = cell_type,
                                 color = consistency
                                 )
                ) +
                scale_color_viridis_c(option = "D") +
                labs(color = "Consistency") +
                geom_point(size = 6.0,
                           alpha = 0.6) +
                geom_label_repel(max.overlaps = Inf,
                                 min.segment.length = 0,
                                 size = 9.0,
                                 force = 6.0) +
                theme(text = element_text(size = 60)) +
                labs(x = cord_x,
                     y = cord_y) +
                theme(legend.key.height = unit(1.5, "cm")) +
                theme(legend.key.width = unit(1.5, "cm"))
            ggsave(filename = paste0(out.dir,
                                     "/figures/consistency_",
                                     con_name,
                                     ".png"),
                   plot = gg_con,
                   dpi = 100,
                   width = 25.0,
                   height = 20.0,
                   limitsize = FALSE)
        }
    }
    ##### Class (plot each label)#####
    if(!is.null(object@eval$external_label)){
        # plot
        label_list <- object@eval$external_label
        for (x in seq_along(label_list)) {
            label_name <- names(label_list[x])
            label_value <- label_list[[x]]
            # convert NA(char) to NA
            label_value[which(label_value=="NA")] <- NA
            df_label <- data.frame(cord_1 = twoD[,1],
                                   cord_2 = twoD[,2],
                                   cell_type = object@union_cellnames,
                                   Class = label_value,
                                   stringsAsFactors = FALSE
                                   )
            gg_label <- ggplot(df_label,
                               aes(x = cord_1,
                                   y = cord_2,
                                   label = cell_type,
                                   color = factor(Class)
                                   )
                             ) +
                labs(color = "Class") +
                geom_point(size = 6.0,
                           alpha = 0.6) +
                geom_label_repel(max.overlaps = Inf,
                                 min.segment.length = 0,
                                 size = 9.0,
                                 force = 6.0) +
                theme(text = element_text(size = 60)) +
                labs(x = cord_x,
                     y = cord_y)
            ggsave(filename = paste0(out.dir,
                                     "/figures/Class_",
                                     label_name,
                                     ".png"),
                   plot = gg_label,
                   dpi = 100,
                   width = 25.0,
                   height = 20.0,
                   limitsize = FALSE
                   )
        }
    }
    ##### 3. 重み/ARIと同定細胞数の関係（例: 論文 Figure 6a）#####
    if(!is.null(object@eval$each_animal)){
        df_each <- object@eval$each_animal
        # sort by weight
        df_sort_weight <- df_each[order(df_each$weight, decreasing=T), ]
        # cowplot
        g1 <- ggplot(df_sort_weight, aes(x = animals,
                                         y= ARI ,
                                         group=1)) +
            geom_line(color = "red", size= 2) +
            scale_x_discrete(limits=df_sort_weight$animals) +
            theme_half_open() +
            theme(text = element_text(size = 36)) +
            theme(axis.title.y=element_text(colour = "red", size = 36)) +
            geom_smooth(method="lm",
                        size =0.5,
                        se = TRUE,
                        alpha = 0.4,
                        color = "red") +
            theme(legend.position = 'none') +
            xlab("Animal Name") +
            theme(axis.text.x= element_text(size = 36,
                                            angle = 45,
                                            hjust = 1)) +
            theme(axis.text.y= element_text(size = 36),
                  plot.background = element_rect(fill = "white",
                                                 color = NA)
                  )

        g2 <- ggplot(df_sort_weight, aes(x = animals,
                                         y= ann_count,
                                         group=1)) +
            geom_line(color = "black", size= 2) +
            scale_x_discrete(limits=df_sort_weight$animals) +
            scale_y_continuous(position = "right") +
            geom_smooth(method="lm",
                        size =0.5,
                        se = TRUE,
                        alpha = 0.4,
                        color = "black") +
            theme_half_open() +
            theme(text = element_text(size = 36)) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            xlab("Animal Name")  +
            ylab("No. of identified cells") +
            theme(axis.text.y= element_text(size = 36))

        aligned_plots_result <- cowplot::align_plots(g1,
                                                     g2,
                                                     align="hv",
                                                     axis="tblr")
        gg_label_ARI <- cowplot::ggdraw(aligned_plots_result[[1]]) +
            cowplot::draw_plot(aligned_plots_result[[2]])
        ggsave(filename = paste0(out.dir, "/figures/weight_ARI_no.png"),
               plot = gg_label_ARI,
               dpi = 100,
               width = 20.0,
               height = 20.0,
               limitsize = FALSE)
    }
    # Output
    object
    }
)

.check_worm_visualize <- function(object,
                                  out.dir,
                                  seed,
                                  tsne.dims,
                                  tsne.perplexity,
                                  tsne.verbose,
                                  tsne.max_iter,
                                  umap.n_neighbors,
                                  umap.n_components,
                                  silhouette.summary){
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
    stopifnot(is.numeric(seed))
    if(seed <= 0 || !.is_integer(seed)){
        stop("Specify seed as a positive integer.")
    }
    stopifnot(is.numeric(tsne.dims))
    if(tsne.dims <= 0 || !.is_integer(tsne.dims)){
        stop("Specify tsne.dims as a positive integer.")
    }
    stopifnot(is.numeric(tsne.perplexity))
    if(tsne.perplexity <= 0 || !.is_integer(tsne.perplexity)){
        stop("Specify tsne.perplexity as a positive integer.")
    }
    stopifnot(is.logical(tsne.verbose))
    stopifnot(is.numeric(tsne.max_iter))
    if(tsne.max_iter <= 0 || !.is_integer(tsne.max_iter)){
        stop("Specify tsne.max_iter as a positive integer.")
    }
    stopifnot(is.numeric(umap.n_neighbors))
    if(umap.n_neighbors <= 0 || !.is_integer(umap.n_neighbors)){
        stop("Specify umap.n_neighbors as a positive integer.")
    }
    stopifnot(is.numeric(umap.n_components))
    if(umap.n_components <= 0 || !.is_integer(umap.n_components)){
        stop("Specify umap.n_components as a positive integer.")
    }
    stopifnot(is.logical(silhouette.summary))
}
