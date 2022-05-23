# ここにroxygen2のコメントを書く
# worm_membership
setGeneric("worm_membership",
    function(object, k){
    standardGeneric("worm_membership")})

# ここにroxygen2のコメントを書く
# worm_clustering
setGeneric("worm_clustering",
    function(object, num.iter=30, thr=1E-10, verbose=FALSE,
        algorithm=c("MCMI", "OINDSCAL", "CSPA")){
    standardGeneric("worm_clustering")})

# ここにroxygen2のコメントを書く
# worm_evaluate
setGeneric("worm_evaluate",
    function(object, labels=NULL){
    standardGeneric("worm_evaluate")})

# ここにroxygen2のコメントを書く
# worm_visualize
setGeneric("worm_visualize",
    function(object, algorithm=c("tSNE", "UMAP"), out.dir=tempdir()){
    standardGeneric("worm_visualize")})
