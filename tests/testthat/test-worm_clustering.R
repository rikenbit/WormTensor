# Pipe Operation
worm_download()$Ds |>
    as_worm_tensor() |>
        worm_membership(k = 6) -> object

expect_equal(length(object@clustering), 0)

######## MCMI #######
worm_clustering(object, verbose = TRUE) -> object_mcmi

expect_equal(object_mcmi@clustering_algorithm, "MCMI")
expect_equal(
    length(object_mcmi@clustering),
    object_mcmi@n_union_cells
)
expect_equal(
    dim(object_mcmi@factor),
    c(object_mcmi@n_union_cells, object_mcmi@k)
)
expect_equal(
    length(object_mcmi@weight),
    object_mcmi@n_animals
)

######## OINDSCAL #######
worm_clustering(object, algorithm = "OINDSCAL", verbose = TRUE) -> object_oindscal

expect_equal(object_oindscal@clustering_algorithm, "OINDSCAL")
expect_equal(
    length(object_oindscal@clustering),
    object_oindscal@n_union_cells
)
expect_equal(
    dim(object_oindscal@factor),
    c(object_oindscal@n_union_cells, object_oindscal@k)
)

######## CSPA #######
worm_clustering(object, algorithm = "CSPA", verbose = TRUE) -> object_cspa

expect_equal(object_cspa@clustering_algorithm, "CSPA")
expect_equal(
    length(object_cspa@clustering),
    object_cspa@n_union_cells
)
expect_equal(
    dim(object_cspa@consensus),
    rep(object_cspa@n_union_cells, 2)
)
