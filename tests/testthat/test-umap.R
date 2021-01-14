test_that("umap works", {
    x <- matrix(rnorm(1e3), nrow = 100)
    expect_is(densmap(x), "matrix")
    Y_init <- prcomp(x)$x[, 1:2]
    expect_is(densmap(x, Y_init = Y_init), "matrix")
    expect_is(umap(x), "matrix")
})
