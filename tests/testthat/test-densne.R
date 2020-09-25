test_that("error conditions", {
    x <- matrix(rnorm(1e4), nrow = 1000)
    x <- x[1:10, ]
    expect_error(densne(x), "perplexity is too large for the number of samples")
    expect_error(
        densne(x, max_iter=-1),
        "number of iterations should be a positive integer"
    )
    expect_error(
        densne(x, perplexity = -1),
        "perplexity should be a positive number"
    )
    expect_error(
        densne(x, exaggeration_factor = -1),
        "exaggeration_factor should be a positive number"
    )
    expect_error(
        densne(x, momentum = -1),
        "momentum should be a positive number"
    )
    expect_error(
        densne(x, eta = -1),
        "eta should be a positive number"
    )
    expect_error(
        densne(x, final_momentum = -1),
        "final_momentum should be a positive number"
    )
    expect_error(
        densne(x, mom_switch_iter = -1),
        "mom_switch_iter should be a positive integer"
    )
    expect_error(
        densne(x, stop_lying_iter = -1),
        "stop_lying_iter should be a positive integer"
    )
    expect_error(
        densne(x, dens_frac = 1.1),
        "dens_frac should lie in \\[0, 1\\]"
    )
    expect_error(
        densne(x, dens_lambda = 1.1),
        "dens_lambda should lie in \\[0, 1\\]"
    )
    expect_error(
        densne(x, theta = 1.1),
        "theta should lie in \\[0, 1\\]"
    )
})

test_that("densne works", {
    x <- matrix(rnorm(1e4), nrow = 1000)
    expect_is(densne(x), "matrix")
    Y_init <- prcomp(x)$x[, 1:2]
    expect_is(densne(x, Y_init = Y_init), "matrix")
})
