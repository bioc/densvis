#' 
#' @examples 
#' x <- matrix(rnorm(1e4), ncol = 1000)
#' d <- densne(x)
#' t <- Rtsne(t(x))
#' @export
densne <- function(
    data,
    dims = 2,
    perplexity = 50,
    theta = 0.5,
    verbose = FALSE,
    max_iter = 1000,
    Y_init = NULL,
    distance_precomputed = FALSE,
    init = FALSE,
    stop_lying_iter = if (is.null(Y_init)) 250L else 0L, 
    mom_switch_iter = if (is.null(Y_init)) 250L else 0L, 
    momentum = 0.5,
    final_momentum = 0.8,
    eta = 200,
    exaggeration_factor = 12,
    dens_frac = 0.3,
    dens_lambda = 0.1,
    final_dens = FALSE,
    num_threads = 1
  ) {

  out <- densne_cpp(
    X = data,
    no_dims = dims,
    perplexity = perplexity,
    theta = theta,
    verbose = verbose,
    max_iter = max_iter,
    distance_precomputed = distance_precomputed,
    Y_in = matrix(),
    init = init,
    stop_lying_iter = stop_lying_iter,
    mom_switch_iter = mom_switch_iter,
    momentum = momentum,
    final_momentum = final_momentum,
    eta = eta,
    exaggeration_factor = exaggeration_factor,
    dens_frac = dens_frac,
    dens_lambda = dens_lambda,
    final_dens = final_dens,
    num_threads = num_threads
  )
  t(out)
}
