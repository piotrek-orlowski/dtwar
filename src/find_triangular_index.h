#include "RcppArmadillo.h"

using namespace Rcpp;

arma::uvec lower_triangular_index(arma::uword matrix_dimension);
  
arma::uvec upper_triangular_index(arma::uword matrix_dimension);