# Pipe Operation
worm_download()$Ds |>
    as_worm_tensor() -> object

# w/ k=3
worm_membership(object, k=3) -> object_k3

expect_equal(is(object_k3@membership_tensor), "Tensor")
expect_equal(object_k3@k, 3)

# w/ k=4
worm_membership(object, k=4) -> object_k4

expect_equal(is(object_k4@membership_tensor), "Tensor")
expect_equal(object_k4@k, 4)

