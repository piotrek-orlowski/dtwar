library(dplyr)

# load data from more-vrp project  ---  2018-03-19 18:35:47  -----
variance_data <- read.csv("../more-vrps/data/covariance-6-fama-french-factors-daily-2018-03-19.csv")

# filter for mkt, smb, hml
variance_data <- variance_data %>%
   filter(of %in% c("mkt_rf","hml","smb"), with %in% c("mkt_rf","hml","smb"))

# set up initial parameters
number_of_series <- 3

# convert to desired matrix form
variance_data_matrix <- variance_data %>% 
  select(date, of, with, cov) %>% 
  arrange(date, with, of) %>% 
  mutate(of_with = paste0(with,"_",of)) %>% 
  select(-of,-with) %>% 
  tidyr::spread(of_with,cov)

eigenvalue_table <- variance_data_matrix %>% 
  group_by(date) %>% 
  do({
    loc_mat <- matrix(0,number_of_series,number_of_series)
    loc_mat[lower.tri(loc_mat, diag = TRUE)] <- as.matrix(select(.,-date))
    loc_mat <- loc_mat + t(loc_mat)
    diag(loc_mat) <- 0.5*diag(loc_mat)
    loc_eigval <- eigen(loc_mat)$values
    loc_eigval <- as.complex(loc_eigval)
    names(loc_eigval) <- paste0("lambda_",1:number_of_series)
    as.data.frame(as.list(loc_eigval))
  })

variance_data_matrix <- variance_data_matrix %>% 
  select(-date) %>% 
  as.matrix()

# variance_data %>% filter(of==with) %>% 
#   arrange(date,of) %>% 
#   group_by(of) %>% 
#   summarise(var_of = var(cov)) %>% 
#   .$var_of

initial_m <- diag(0.99,number_of_series,number_of_series)

initial_sigma <- 2*diag(variance_data %>% filter(of==with) %>% 
                        arrange(date,of) %>% 
                        group_by(of) %>% 
                        summarise(var_of = var(cov)) %>% 
                        .$var_of)

# evaluate objective
obj_test <- estimating_function_first_order_mm(parameters = c(as.numeric(initial_m), as.numeric(initial_sigma[lower.tri(initial_sigma, diag = TRUE)]))
                                               , data = variance_data_matrix
                                               , model_dimension = number_of_series)


# try optimisation
library(nloptr)

# set some parameter bounds
sigma_lower_bound_mat <- initial_sigma
diag(sigma_lower_bound_mat) <- 0
sigma_lower_bound_mat[lower.tri(sigma_lower_bound_mat)] <- -0.3
sigma_upper_bound_mat <- sigma_lower_bound_mat * (-1)
diag(sigma_upper_bound_mat) <- 0.5

lower_bounds <- c(rep(-4,number_of_series^2), as.numeric(sigma_lower_bound_mat[lower.tri(sigma_lower_bound_mat, diag = TRUE)]))
upper_bounds <- c(rep(4,number_of_series^2), as.numeric(sigma_upper_bound_mat[lower.tri(sigma_upper_bound_mat, diag = TRUE)]))

optim_res <- nloptr::lbfgs(x0 = c(as.numeric(initial_m), as.numeric(initial_sigma[lower.tri(initial_sigma, diag = TRUE)]))
                           , fn = estimating_function_first_order_mm
                           , nl.info = TRUE
                           , lower = lower_bounds
                           , upper = upper_bounds
                           , control = list(xtol_rel = 1e-12, maxeval = 1e6)
                           , data = variance_data_matrix
                           , model_dimension = number_of_series)

m_estimated <- matrix(optim_res$par[1:number_of_series^2],number_of_series,number_of_series)
sigma_estimated <- matrix(0,number_of_series,number_of_series)
sigma_estimated[lower.tri(sigma_estimated, diag = TRUE)] <- optim_res$par[-c(1:number_of_series^2)]
sigma_estimated <- sigma_estimated + t(sigma_estimated)
diag(sigma_estimated) <- 0.5 * diag(sigma_estimated)

eigen(m_estimated)$values
eigen(sigma_estimated)$values

estimating_function_second_order_mm(portfolio_allocation = c(1,1,1), sigma_star_matrix = sigma_estimated, m_matrix = m_estimated, data = variance_data_matrix)[1]

# try with constraints on eigenvalues
library(alabama)
equality_constraint <- function(parameters, data, model_dimension){
  
  m_matrix <- matrix(parameters[1:model_dimension^2], model_dimension, model_dimension)
  constr_m <- Im(eigen(m_matrix)$values)
  
  s_matrix <- matrix(0, model_dimension, model_dimension)
  s_matrix[lower.tri(s_matrix, diag=TRUE)] <- parameters[-(1:model_dimension^2)]
  s_matrix <- s_matrix + t(s_matrix)
  diag(s_matrix) <- 0.5 * diag(s_matrix)
  
  constr_s <- Im(eigen(s_matrix)$values)
  
  return(c(constr_m, constr_s))
}

inequality_constraint <- function(parameters, data, model_dimension){
  
  m_matrix <- matrix(parameters[1:model_dimension^2], model_dimension, model_dimension)
  constr_m <- Re(eigen(m_matrix)$values)
  
  s_matrix <- matrix(0, model_dimension, model_dimension)
  s_matrix[lower.tri(s_matrix, diag=TRUE)] <- parameters[-(1:model_dimension^2)]
  s_matrix <- s_matrix + t(s_matrix)
  diag(s_matrix) <- 0.5 * diag(s_matrix)
  
  constr_s <- Re(eigen(s_matrix)$values)
  
  return(c(constr_m, constr_s))
}

alabama_res <- alabama::auglag(par = c(as.numeric(initial_m), as.numeric(initial_sigma[lower.tri(initial_sigma, diag = TRUE)]))
                      , fn = estimating_function_first_order_mm
                      , hin = inequality_constraint
                      , heq = equality_constraint
                      , control.outer = list(method = "nlminb", kkt2.check=FALSE, eps = 1e-9)
                      , control.optim = list(rel.tol = 1e-12, abs.tol = 1e-20, x.tol = 1e-12)
                      , data = variance_data_matrix
                      , model_dimension = number_of_series
                      )


# try again with constraints on eigenvalues
library(nloptr)
equality_constraint <- function(parameters){
  
  model_dimension <- number_of_series
  
  m_matrix <- matrix(parameters[1:model_dimension^2], model_dimension, model_dimension)
  constr_m <- Im(eigen(m_matrix)$values)
  
  s_matrix <- matrix(0, model_dimension, model_dimension)
  s_matrix[lower.tri(s_matrix, diag=TRUE)] <- parameters[-(1:model_dimension^2)]
  s_matrix <- s_matrix + t(s_matrix)
  diag(s_matrix) <- 0.5 * diag(s_matrix)
  
  constr_s <- Im(eigen(s_matrix)$values)
  
  return(c(constr_m, constr_s))
}

inequality_constraint <- function(parameters){
  
  model_dimension <- number_of_series
  
  m_matrix <- matrix(parameters[1:model_dimension^2], model_dimension, model_dimension)
  constr_m <- Re(eigen(m_matrix)$values)
  
  return(c(constr_m))
}

auglag_res <- nloptr::auglag(x0 = c(as.numeric(initial_m), as.numeric(initial_sigma[lower.tri(initial_sigma, diag = TRUE)]))
                      , fn = estimating_function_first_order_mm
                      , hin = inequality_constraint
                      , heq = equality_constraint
                      , localsolver = "LBFGS"
                      , control = list(xtol_rel = 1e-12, maxeval = 1e6, ftol_abs = 1e-20, ftol_rel = 1e-18)
                      , data = variance_data_matrix
                      , model_dimension = number_of_series
)

m_estimated_auglag <- matrix(auglag_res$par[1:number_of_series^2],number_of_series,number_of_series)
sigma_estimated_auglag <- matrix(0,number_of_series,number_of_series)
sigma_estimated_auglag[lower.tri(sigma_estimated, diag = TRUE)] <- auglag_res$par[-c(1:number_of_series^2)]
sigma_estimated_auglag <- sigma_estimated_auglag + t(sigma_estimated_auglag)
diag(sigma_estimated_auglag) <- 0.5 * diag(sigma_estimated_auglag)

eigen(m_estimated_auglag)$values
eigen(sigma_estimated_auglag)$values


estimating_function_second_order_mm(portfolio_allocation = c(1,1,1), sigma_star_matrix = sigma_estimated_auglag, m_matrix = m_estimated_auglag, data = variance_data_matrix)[1]

