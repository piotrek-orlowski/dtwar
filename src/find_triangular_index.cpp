#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::uvec lower_triangular_index(arma::uword matrix_dimension){
  
  // initialize matrix
  arma::mat square_matrix(matrix_dimension, matrix_dimension, arma::fill::ones);
  
  // lower-triangular part of sigma-star
  arma::uvec triangular_part = arma::find(arma::trimatl(square_matrix));
  
  return triangular_part;
}

arma::uvec upper_triangular_index(arma::uword matrix_dimension){
  
  // initialize matrix
  arma::mat square_matrix(matrix_dimension, matrix_dimension, arma::fill::ones);
  
  // lower-triangular part of sigma-star
  arma::uvec triangular_part = arma::find(arma::trimatu(square_matrix));
  
  return triangular_part;
}