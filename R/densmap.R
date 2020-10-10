#' Density-preserving t-SNE
#' @param x A numeric matrix or matrix-like object.
#' @param n_components The dimension of the space to embed 
#' into. This defaults to 2 to provide easy visualization, 
#' but can reasonably be set to any integer value in the 
#' range 2 to 100.
#' @param dens_frac numeric; fraction of the iterations for
#' which the full objective function (including the 
#' density-preserving term) is used. For the first 
#' \code{1 - dens_frac} fraction of the iterations, only
#' the original t-SNE objective function is used.
#' @param dens_lambda numeric; the relative importance of the 
#' density-preservation term compared to the original t-SNE 
#' objective function.
#' @param var_shift Regularization term added to the variance 
#' of embedding local radius for stability (float, 
#' non-negative); default 0.1.
#' @param n_neighbors The size of local neighborhood 
#' (in terms of number of neighboring sample points) used for 
#' manifold approximation. Larger values result in more 
#' global views of the manifold, while smaller values result 
#' in more local data being preserved. In general values 
#' should be in the range 2 to 100.
#' @param metric The metric to use to compute distances in 
#' high dimensional space. If a string is passed it must match
#' one of:
#' - "euclidean"
#' - "manhattan"
#' - "chebyshev"
#' - "minkowski"
#' - "canberra"
#' - "braycurtis"
#' - "mahalanobis"
#' - "wminkowski"
#' - "seuclidean"
#' - "cosine"
#' - "correlation"
#' - "haversine"
#' - "hamming"
#' - "jaccard"
#' - "dice"
#' - "russelrao"
#' - "kulsinski"
#' - "rogerstanimoto"
#' - "sokalmichener"
#' - "sokalsneath"
#' - "yule"
#' @param n_epochs The number of training epochs to be used 
#' in optimizing the low dimensional embedding. Larger values 
#' result in more accurate embeddings. If None is specified a 
#' value will be selected based on the size of the input 
#' dataset (200 for large datasets, 500 for small).
#' a valid predefined metric.
#' @param learning_rate The initial learning rate for the embedding optimization.
#' @param init How to initialize the low dimensional 
#' embedding. Valid options:
#' - "spectral": use a spectral embedding of the fuzzy 1-skeleton
#' - "random": assign initial embedding positions at random.
#' @param min_dist The effective minimum distance between 
#' embedded points. Smaller values will result in a more 
#' clustered/clumped embedding where nearby points on the 
#' manifold are drawn closer together, while larger values 
#' will result on a more even dispersal of points. The value 
#' should be set relative to the spread value, which 
#' determines the scale at which embedded points will be 
#' spread out.
#' @param spread The effective scale of embedded points. In combination with min_dist this determines how clustered/clumped the embedded points are.
#' @param set_op_mix_ratio Interpolate between (fuzzy) union 
#' and intersection as the set operation used to combine 
#' local fuzzy simplicial sets to obtain a global fuzzy 
#' simplicial sets. Both fuzzy set operations use the product 
#' t-norm. The value of this parameter should be between 0.0 
#' and 1.0; a value of 1.0 will use a pure fuzzy union, while 
#' 0.0 will use a pure fuzzy intersection.
#' @param local_connectivity The local connectivity required 
#' – i.e. the number of nearest neighbors that should be 
#' assumed to be connected at a local level. The higher this 
#' value the more connected the manifold becomes locally. In 
#' practice this should be not more than the local intrinsic 
#' dimension of the manifold.
#' @param repulsion_strength Weighting applied to negative 
#' samples in low dimensional embedding optimization. Values 
#' higher than one will result in greater weight being given 
#' to negative samples.
#' @param negative_sample_rate The number of negative samples 
#' to select per positive sample in the optimization process. 
#' Increasing this value will result in greater repulsive 
#' force being applied, greater optimization cost, but 
#' slightly more accuracy.
#' @param transform_queue_size For transform operations 
#' (embedding new points using a trained model_ this will 
#' control how aggressively to search for nearest neighbors. 
#' Larger values will result in slower performance but more 
#' accurate nearest neighbor evaluation.
#' @param random_state The seed used by the random number 
#' generator.
#' @param angular_rp_forest Whether to use an angular random 
#' projection forest to initialise the approximate nearest 
#' neighbor search. This can be faster, but is mostly on 
#' useful for metric that use an angular style distance such 
#' as cosine, correlation etc. In the case of those metrics 
#' angular forests will be chosen automatically.
#' @param target_n_neighbors The number of nearest neighbors 
#' to use to construct the target simplcial set. If set to -1 
#' use the n_neighbors value.
#' @param target_weight Weighting factor between data 
#' topology and target topology. A value of 0.0 weights 
#' entirely on data, a value of 1.0 weights entirely on 
#' target. The default of 0.5 balances the weighting equally 
#' between data and target.
#' @return A numeric matrix
#' @examples
#' set.seed(42)
#' x <- matrix(rnorm(2000), ncol=20) 
#' densmap(x)
#' @export
densmap <- function(
        x,
        n_components = 2L,
        dens_frac = 0.3,
        dens_lambda = 0.1,
        var_shift = 0.1,
        n_neighbors = 15L,
        metric = "euclidean",
        n_epochs = 750L,
        learning_rate = 1.0,
        init = c("spectral", "random"),
        min_dist = 0.1,
        spread = 0.1,
        # low_memory = FALSE,
        set_op_mix_ratio = 1.0,
        local_connectivity = 1L,
        repulsion_strength = 1.0,
        negative_sample_rate = 5L,
        transform_queue_size = 4.0,
        random_state = NULL,
        angular_rp_forest = FALSE,
        target_n_neighbors = -1,
        target_weight = 0.5
    ) {
    x <- as.matrix(x)
    init <- match.arg(init)
    .checks(
        x = x,
        n_components = n_components,
        dens_frac = dens_frac,
        dens_lambda = dens_lambda,
        n_neighbors = n_neighbors,
        metric = metric,
        n_epochs = n_epochs,
        var_shift = var_shift,
        learning_rate = learning_rate,
        init = init,
        min_dist = min_dist,
        spread = spread,
        # low_memory = low_memory,
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        repulsion_strength = repulsion_strength,
        negative_sample_rate = negative_sample_rate,
        transform_queue_size = transform_queue_size,
        random_state = random_state,
        angular_rp_forest = angular_rp_forest,
        target_n_neighbors = target_n_neighbors,
        target_weight = target_weight
    )
    proc <- basiliskStart(python_env)
    on.exit(basiliskStop(proc))
    out <- basiliskRun(proc, .fit_densmap,
        x = x, 
        n_components = n_components,
        dens_frac = dens_frac,
        dens_lambda = dens_lambda,
        n_neighbors = n_neighbors,
        metric = metric,
        n_epochs = n_epochs,
        var_shift = var_shift,
        learning_rate = learning_rate,
        init = init,
        min_dist = min_dist,
        spread = spread,
        # low_memory = low_memory,
        set_op_mix_ratio = set_op_mix_ratio,
        local_connectivity = local_connectivity,
        repulsion_strength = repulsion_strength,
        negative_sample_rate = negative_sample_rate,
        transform_queue_size = transform_queue_size,
        random_state = random_state,
        angular_rp_forest = angular_rp_forest,
        target_n_neighbors = target_n_neighbors,
        target_weight = target_weight

    )
    out

}

.fit_densmap <- function(x, ...) {
    densmap <- reticulate::import("densmap")
    densmap$densMAP(
        final_dens = FALSE,
        ...
    )$fit_transform(x)
}


.checks <- function(
        x,
        n_components,
        dens_frac,
        dens_lambda,
        n_neighbors,
        metric,
        n_epochs,
        var_shift,
        learning_rate,
        init,
        min_dist,
        spread,
        low_memory,
        set_op_mix_ratio,
        local_connectivity,
        repulsion_strength,
        negative_sample_rate,
        transform_queue_size,
        random_state,
        angular_rp_forest,
        target_n_neighbors,
        target_weight
    ) {
    match.arg(metric,
        choices = c(
            "euclidean",
            "manhattan",
            "chebyshev",
            "minkowski",
            "canberra",
            "braycurtis",
            "mahalanobis",
            "wminkowski",
            "seuclidean",
            "cosine",
            "correlation",
            "haversine",
            "hamming",
            "jaccard",
            "dice",
            "russelrao",
            "kulsinski",
            "rogerstanimoto",
            "sokalmichener",
            "sokalsneath",
            "yule"
        )
    )
    assert_that(
        is.integer(n_components),
        n_components > 0,
        length(n_components) == 1,
        is.numeric(dens_frac),
        dens_frac >= 0,
        dens_frac <= 1,
        is.numeric(dens_lambda),
        dens_lambda >= 0,
        dens_lambda <= 1,
        n_neighbors > 0,
        is.integer(n_neighbors),
        length(n_neighbors) == 1,
        is.integer(n_epochs),
        n_epochs > 0,
        length(n_epochs) == 1,
        is.numeric(var_shift),
        var_shift > 0,
        length(var_shift) == 1,
        is.numeric(learning_rate),
        learning_rate > 0,
        length(learning_rate) == 1,
        is.numeric(min_dist),
        min_dist > 0,
        length(min_dist) == 1,
        is.numeric(spread),
        spread > 0,
        length(spread) == 1,
        is.numeric(set_op_mix_ratio),
        set_op_mix_ratio > 0,
        length(set_op_mix_ratio) == 1,
        is.integer(local_connectivity),
        local_connectivity > 0,
        length(local_connectivity) == 1,
        is.numeric(repulsion_strength),
        repulsion_strength > 0,
        length(repulsion_strength) == 1,
        is.integer(negative_sample_rate),
        negative_sample_rate > 0,
        length(negative_sample_rate) == 1,
        is.numeric(transform_queue_size),
        transform_queue_size > 0,
        length(transform_queue_size) == 1,
        is.null(random_state) ||
            is.integer(random_state)
            & random_state > 0
            & length(random_state) == 1,
        is.logical(angular_rp_forest),
        length(angular_rp_forest) == 1,
        is.numeric(target_n_neighbors),
        is.numeric(target_weight),
        length(target_weight) == 1
    )
}