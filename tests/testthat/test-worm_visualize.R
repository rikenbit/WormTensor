# Temporary directory to save figures
out.dir <- tempdir()

# Pipe Operation
worm_download()$Ds |>
    as_worm_tensor() |>
        worm_membership(k=3) |>
            worm_clustering() |>
                worm_evaluate() |>
                    worm_visualize(out.dir) -> object

# ここにテストを書いていく
filename1 <- paste0(out.dir, "/XXX")
expect_true(file.exists(filename1))

filename2 <- paste0(out.dir, "/XXX")
expect_true(file.exists(filename2))

filename3 <- paste0(out.dir, "/XXX")
expect_true(file.exists(filename3))

filename4 <- paste0(out.dir, "/XXX")
expect_true(file.exists(filename4))

filename5 <- paste0(out.dir, "/XXX")
expect_true(file.exists(filename5))
