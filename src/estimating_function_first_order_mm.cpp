#include "RcppArmadillo.h"
#include <stdexcept>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double estimating_function_first_order_mm(const arma::colvec & parameters
                                            , const arma::mat & data
                                            , const unsigned int model_dimension
) {
  
  // initialize estimation criterion
  double objective = 0.0;
  
  // parameter vector lengths
  unsigned int length_for_m = model_dimension * model_dimension;
  unsigned int length_for_sigma = model_dimension * (model_dimension+1) / 2;
  
  // check model dimensions
  unsigned int required_parameter_length = length_for_m + length_for_sigma;
  
  if(required_parameter_length != parameters.n_elem){
    throw std::invalid_argument("Incorrect length of parameter vector");
  }
  
  if(data.n_cols != length_for_sigma){
    throw std::invalid_argument("Incorrect number of columns of data matrix");
  }
  
  // initialize m_matrix and reshape
  arma::mat m_matrix;
  m_matrix.insert_cols(0, parameters.subvec(0,length_for_m - 1));
  m_matrix.reshape(model_dimension, model_dimension);
  
  // initialize sigma_star, fill with ones
  arma::mat sigma_star_matrix(model_dimension, model_dimension, arma::fill::ones);
  
  // lower-triangular part of sigma-star
  arma::uvec triangular_part_of_sigma = arma::find(arma::trimatl(sigma_star_matrix));
  
  sigma_star_matrix(triangular_part_of_sigma) = parameters.subvec(length_for_m, length_for_m + length_for_sigma - 1);
  sigma_star_matrix = arma::symmatl(sigma_star_matrix);
  
  // initialize variables for loop
  arma::mat y_matrix_t(model_dimension, model_dimension, arma::fill::zeros);
  arma::mat y_matrix_t_minus(model_dimension, model_dimension, arma::fill::zeros);
  
  // determine number of obserations
  unsigned int number_of_obs = data.n_rows;
  
  // loop over data to calculate criterion
  // start at 2nd (=1) row
  for(unsigned int row_number = 1; row_number < number_of_obs; ++row_number){
    
    y_matrix_t(triangular_part_of_sigma) = data.row(row_number).t();
    y_matrix_t = arma::symmatl(y_matrix_t);
    
    // if(row_number == 1){
    //   y_matrix_t.print("y_2");
    // }
    
    y_matrix_t_minus(triangular_part_of_sigma) = data.row(row_number-1).t();
    y_matrix_t_minus = arma::symmatl(y_matrix_t_minus);
    
    // if(row_number == 1){
    //   y_matrix_t_minus.print("y_1");
    // }
    
    objective += pow(arma::norm(arma::vectorise(y_matrix_t) - arma::vectorise(m_matrix * y_matrix_t_minus * m_matrix.t() + sigma_star_matrix)),2.0);
  }
  
  return objective;
}