# Temporary directory to save figures
out.dir <- tempdir()

# Pipe Operation
worm_download()$Ds |>
    as_worm_tensor() |>
        worm_membership(k=6) |>
            worm_clustering() |>
                worm_evaluate() |>
                    worm_visualize() -> object

######### no labels #########
filename1 <- paste0(out.dir, "/figures/Silhouette.png")
expect_true(file.exists(filename1))

filename2 <- paste0(out.dir, "/figures/Cluster.png")
expect_true(file.exists(filename2))

filename3 <- paste0(out.dir, "/figures/no_identified.png")
expect_true(file.exists(filename3))
