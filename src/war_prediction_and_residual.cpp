#include "RcppArmadillo.h"
#include "find_triangular_index.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List war_prediction_and_residuals(const arma::mat m_matrix, const arma::mat sigma_matrix, const arma::mat data){

  // initialize matrices  
  arma::mat prediction_matrix(data);
  arma::mat residual_matrix(data);
  
  prediction_matrix.zeros();
  residual_matrix.zeros();
  
  arma::mat y_temp_matrix(m_matrix.n_rows, m_matrix.n_cols, arma::fill::zeros);
  
  arma::uvec lower_triangular = lower_triangular_index(y_temp_matrix.n_rows);
  
  y_temp_matrix(lower_triangular) = data.row(0).t();
  y_temp_matrix = symmatl(y_temp_matrix);
  
  arma::mat prediction_temp_matrix(y_temp_matrix);
  prediction_temp_matrix.zeros();
  
  // determine number of obserations
  unsigned int number_of_obs = data.n_rows;
  
  // calculate predictions and residuals
  for(unsigned int row_number = 1; row_number < number_of_obs; row_number++){
    
    // predict
    prediction_temp_matrix = m_matrix * y_temp_matrix * m_matrix.t() + sigma_matrix;
    
    // save prediction
    prediction_matrix.row(row_number) = arma::vectorise(prediction_temp_matrix(lower_triangular)).t();
    
    // calculate and save residuals
    residual_matrix.row(row_number) = data.row(row_number) - prediction_matrix(row_number);
    
    // write next obs to y_temp
    y_temp_matrix(lower_triangular) = data.row(row_number).t();
    y_temp_matrix = symmatl(y_temp_matrix);
  }
  
  List result_list;
  result_list = List::create(Rcpp::Named("predicted") = prediction_matrix
                               , Rcpp::Named("residuals") = residual_matrix
                               );
  
  return result_list;
}