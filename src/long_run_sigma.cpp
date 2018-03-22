#include "RcppArmadillo.h"
#include <stdexcept>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat long_run_sigma(arma::mat m_matrix, arma::mat sigma_star_matrix){
  
  // infer dimension
  unsigned int model_dimension = m_matrix.n_rows;
  
  // check dim compatibility
  if(model_dimension != m_matrix.n_cols){
    throw std::invalid_argument("m_matrix not square");
  }
  
  if((model_dimension != sigma_star_matrix.n_cols) | (model_dimension != sigma_star_matrix.n_rows)){
    throw std::invalid_argument("wrong sigma_star dimension");
  }
  
  // attempt analytical solution
  arma::mat unit_n_squared(model_dimension * model_dimension, model_dimension * model_dimension, arma::fill::eye);
  arma::mat m_kron_m = arma::kron(m_matrix, m_matrix);
  arma::vec sigma_vec = arma::vectorise(sigma_star_matrix);

  arma::mat sigma_long_vec = arma::solve(unit_n_squared - m_kron_m,sigma_vec);

  sigma_long_vec.reshape(model_dimension, model_dimension);
  
  return sigma_long_vec;
}