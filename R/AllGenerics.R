# ここにroxygen2のコメントを書く
# worm_membership
setGeneric("worm_membership",
    function(object, k){
    standardGeneric("worm_membership")})

# ここにroxygen2のコメントを書く
# worm_clustering
setGeneric("worm_clustering",
    function(object,
             num.iter=30,
             thr=1E-10,
             verbose=FALSE,
             algorithm=c("MCMI", "OINDSCAL", "CSPA")
             ){
    standardGeneric("worm_clustering")})

# ここにroxygen2のコメントを書く
# worm_evaluate
setGeneric("worm_evaluate",
    function(object, labels=NULL){
    standardGeneric("worm_evaluate")})

# ここにroxygen2のコメントを書く
# worm_visualize
setGeneric("worm_visualize",
    function(object,
             algorithm=c("tSNE", "UMAP"),
             out.dir=tempdir(),
             seed=1234,
             tsne.dims=2,
             tsne.perplexity=15,
             tsne.verbose=FALSE,
             tsne.max_iter=1000,
             umap.n_neighbors=15,
             umap.n_components=2,
             silhouette.summary=FALSE){
    standardGeneric("worm_visualize")})
