#include <Rcpp.h>
#include "densne.h"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix densne_cpp(
    NumericMatrix X,
    int no_dims,
    double perplexity,
    double theta,
    bool verbose,
    int max_iter,
    NumericMatrix Y_in,
    bool init,
    int stop_lying_iter,
    int mom_switch_iter,
    double momentum,
    double final_momentum,
    double eta,
    double exaggeration_factor,
    double dens_frac,
    double dens_lambda,
    bool final_dens,
    unsigned int num_threads
  ) {


  size_t N = X.ncol(), D = X.nrow();
  double * data=X.begin();
  
  if (verbose) {
    Rprintf("Read the %i x %i data matrix successfully!\n", N, D);
  }

  std::vector<double> Y(N * no_dims), costs(N), itercosts(static_cast<int>(std::ceil(max_iter/50.0)));

  // Providing user-supplied solution.
  if (init) {
      for (size_t i = 0; i < Y.size(); i++) Y[i] = Y_in[i];
      if (verbose) Rprintf("Using user supplied starting positions\n");
  }
  
  double* dens = NULL;
  if (final_dens) {
    dens = (double*) malloc(N * 2 * sizeof(double));
  }

  // Run densne
  DENSNE::run(
    data,
    N,
    D,
    Y.data(),
    dens,
    no_dims,
    perplexity,
    theta,
    false,
    max_iter,
    momentum,
    final_momentum,
    eta,
    stop_lying_iter,
    mom_switch_iter,
    dens_frac,
    dens_lambda,
    final_dens,
    verbose
  );

  // double* X,
  // int N,
  // int D,
  // double* Y,
  // double* dens,
  // int no_dims,
  // double perplexity,
  // double theta, 
  // bool skip_random_init,
  // int max_iter,
  // int stop_lying_iter,
  // int mom_switch_iter,
  // double dens_frac,
  // double dens_lambda,
  // bool final_dens

  return Rcpp::NumericMatrix(no_dims, N, Y.data());
  // todo: make a matrix from dens
  // Rcpp::List(
  //   Rcpp::Named("x") = Rcpp::NumericMatrix(no_dims, N, Y.data()),
  //   Rcpp::Named("dens") = Rcpp::NumericMatrix(no_dims, N, dens)
  // );

  // return Rcpp::List::create(Rcpp::_["Y"]=Rcpp::NumericMatrix(no_dims, N, Y.data()), 
  //         Rcpp::_["costs"]=Rcpp::NumericVector(costs.begin(), costs.end()),
  //         Rcpp::_["itercosts"]=Rcpp::NumericVector(itercosts.begin(), itercosts.end()));
}
