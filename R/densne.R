#' Density-preserving t-SNE
#' @param data Input data matrix.
#' @param dims Integer output dimensionality.
#' @param perplexity Perplexity parameter (should not be bigger than 
#'  3 * perplexity < nrow(X) - 1).
#' @param theta Speed/accuracy trade-off (increase for less
#'  accuracy), set to 0.0 for exact TSNE 
#' @param verbose Logical; Whether progress updates should be printed 
#' @param max_iter integer; Number of iterations
#' @param Y_init matrix; Initial locations of the objects. If NULL, random
#' initialization will be used 
#' @param stop_lying_iter integer; Iteration after which the perplexities are no
#'  longer exaggerated 
#' @param mom_switch_iter integer; Iteration after which the final momentum is
#'  used
#' @param momentum numeric; Momentum used in the first part of the optimization
#' @param final_momentum numeric; Momentum used in the final part of the 
#'  optimization
#' @param eta numeric; Learning rate
#' @param exaggeration_factor numeric; Exaggeration factor used to multiply the
#'  affinities matrix P in the first part of the optimization
#' @param dens_frac numeric; fraction of the iterations for which the full 
#'  objective function (including the density-preserving term) is used.
#'  For the first \code{1 - dens_frac} fraction of the iterations, only
#'  the original t-SNE objective function is used.
#' @param dens_lambda numeric; the relative importanceof the 
#'  density-preservation term compared to the original t-SNE objective function.
#' @param num_threads Number of threads to be used for parallelisation.
#' @return A numeric matrix corresponding to the t-SNE embedding 
#' @references
#' Density-Preserving Data Visualization Unveils Dynamic Patterns of Single-Cell 
#' Transcriptomic Variability
#' Ashwin Narayan, Bonnie Berger, Hyunghoon Cho;
#' bioRxiv (2020)
#' <doi:10.1101/2020.05.12.077776>
#' @examples 
#' x <- matrix(rnorm(1e4), nrow = 1000)
#' d <- densne(x)
#' plot(d)
#' @export
densne <- function(
        data,
        dims = 2,
        perplexity = 50,
        theta = 0.5,
        verbose = getOption("verbose", FALSE),
        max_iter = 1000,
        Y_init = NULL,
        stop_lying_iter = if (is.null(Y_init)) 250L else 0L, 
        mom_switch_iter = if (is.null(Y_init)) 250L else 0L, 
        momentum = 0.5,
        final_momentum = 0.8,
        eta = 200,
        exaggeration_factor = 12,
        dens_frac = 0.3,
        dens_lambda = 0.1,
        num_threads = 1
    ) {
    data <- t(as.matrix(data))

    out <- densne_cpp(
        X = data,
        no_dims = dims,
        perplexity = perplexity,
        theta = theta,
        verbose = verbose,
        max_iter = max_iter,
        Y_in = matrix(),
        init = !is.null(Y_init),
        stop_lying_iter = stop_lying_iter,
        mom_switch_iter = mom_switch_iter,
        momentum = momentum,
        final_momentum = final_momentum,
        eta = eta,
        exaggeration_factor = exaggeration_factor,
        dens_frac = dens_frac,
        dens_lambda = dens_lambda,
        final_dens = FALSE,
        num_threads = num_threads
    )
    t(out)
}
