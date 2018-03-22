#include "RcppArmadillo.h"
#include "find_triangular_index.h"
#include "long_run_sigma.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec estimating_function_second_order_mm(const arma::colvec portfolio_allocation
                                                , const arma::mat sigma_star_matrix
                                                , const arma::mat m_matrix
                                                , const arma::mat data
){
  
  unsigned int model_dimension = m_matrix.n_rows;
    
  arma::mat long_run_sigma_matrix;
  
  long_run_sigma_matrix = long_run_sigma(m_matrix, sigma_star_matrix);
  
  arma::mat y_temp_matrix(m_matrix.n_rows, m_matrix.n_cols, arma::fill::zeros);
  
  arma::uvec lower_triangular = lower_triangular_index(y_temp_matrix.n_rows);
  
  // total volatility
  arma::vec portfolio_volatilities(data.n_rows, arma::fill::zeros);
  double volatility_estimate;
  
  // determine number of obserations
  unsigned int number_of_obs = data.n_rows;
  
  // calculate predictions and residuals
  for(unsigned int row_number = 0; row_number < number_of_obs; row_number++){
    
    y_temp_matrix(lower_triangular) = data.row(row_number).t();
    y_temp_matrix = arma::symmatl(y_temp_matrix);
      
    portfolio_volatilities(row_number) = arma::as_scalar(portfolio_allocation.t() * y_temp_matrix * portfolio_allocation);
  }
  
  portfolio_volatilities = portfolio_volatilities - arma::as_scalar(arma::mean(portfolio_volatilities));
  portfolio_volatilities = arma::pow(portfolio_volatilities,2.0);
  
  volatility_estimate = arma::as_scalar(arma::mean(portfolio_volatilities));
  
  double k_estimate = 2.0 * arma::as_scalar(arma::pow(portfolio_allocation.t() * long_run_sigma_matrix * portfolio_allocation, 2.0)) / volatility_estimate;
  
  arma::mat sigma_estimate = 1.0/k_estimate * sigma_star_matrix;
  
  arma::vec output_vector(1 + model_dimension*(model_dimension+1)/2,arma::fill::zeros);
  
  output_vector(0) = k_estimate;
  output_vector.tail(model_dimension*(model_dimension+1)/2) = arma::vectorise(sigma_estimate(lower_triangular));
  
  return output_vector;
}