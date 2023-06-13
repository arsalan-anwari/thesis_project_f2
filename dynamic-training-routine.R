
# Load libraries 

library(sf) # used to load maps and manipulate spatial features
library(tidyverse) # used for piping functionality of data manipulation functions
library(gstat) # used for interpolation and validation of idw and kriging models
library(caret) # used for cross validation of other models

# Set working directory
setwd("C:/Users/arsal/Workspaces/thesis_project_f2")

# Set random seed
set.seed(73829)

# Load data
sf_observations <- readRDS("data/dynamic-training-routine.rds")

# Split x and y coordinate from observations and receptors as separate columns
sf_observations <- sf_observations %>%
  mutate(
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  ) %>%
  select(x, y, windspeed)

# Define the function used to calculate training weights based on RSME gain/loss
loss_function <- function(rsme_old, rsme_new, power=1.09748732){
  rsme_diff <- abs(rsme_old - rsme_new)
  rsme_loss <- -rsme_diff^power
  return( abs(rsme_loss) )
}

# Define the function used to dynamically train a single variation of a model
dynamic_model_train <- function(rsme_calc_func, model_varient, base_weight = 0.0, loss_threshold = 0.0005, max_iterations = 20){
  
  # Pre-calculate the first two RSME values.
  old_rsme <- rsme_calc_func(0, model_varient)
  new_rsme <- rsme_calc_func(loss_threshold, model_varient)
  
  # Pre-calculate the first loss_value
  loss_value <- loss_function(old_rsme, new_rsme)
  
  # Save states of pre-calculated values
  best_fit_weight <- loss_value
  old_rsme <- new_rsme
  
  for( i in 0:max_iterations ){
    new_rsme <- rsme_calc_func(best_fit_weight, model_varient)
    loss_value <- loss_function(old_rsme, new_rsme)
    best_fit_weight <- best_fit_weight + loss_value
    
    if(loss_value < loss_threshold){ return(base_weight + best_fit_weight) }
    
    old_rsme <- new_rsme
  }
  
  return(base_weight + best_fit_weight)
  
}


# Define the function used to dynamically train a all variation of a model
dynamic_model_train_vars <- function(rsme_calc_func, model_variants, base_weight = 0.0, loss_threshold = 0.0005, max_iterations = 20){
  n <- length(model_variants)
  best_fit_weights <- rep(0.01, n)
  
  for (i in 1:n) {
    model_variant <- model_variants[i]
    best_fit_weights[i] <- dynamic_model_train(rsme_calc_func, model_variant, base_weight, loss_threshold, max_iterations)
  }
  
  return(best_fit_weights)
}

# Define the function used to calculate the RSME from residuals
rmse <- function(residuals) {
  return(sqrt(sum((residuals)^2)/length(residuals)))
}

# Define the functions and helper functions used to calculate the RSME for each model

## Trend surface
trend_surface.calc_new_rsme <- function(weight, formula = "windspeed ~ x + y"){
  model_cv <- train(
    as.formula(formula),
    method = "lm",
    data = sf_observations,
    trControl = trainControl(method = "LOOCV"),
    weights = (1 / windspeed) + weight
  )
  
  return( model_cv$results$RMSE )
}

## MQ-RBF
mqrbf.base_weight <- 1.0

mqrbf.calc_alpha_weights <- function(x, y, z, eps) {
  n <- length(x)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- sqrt(1 + eps * ((sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2))^2) ) 
    }
  }
  alpha <- solve(K) %*% z
  return(as.vector(alpha))
}

mqrbf.calc_distance <- function(r, c){
  return ( sqrt( (r^2) + (c^2) ) )
}

mqrbf.predict <- function(
    x_receptors, y_receptors, 
    x_observations, y_observations,
    alpha_weights, smoothing_factor
){
  
  n_i <- length(x_receptors)
  n_j <- length(x_observations)
  
  # value of 0.05 is to ensure that a float is used and not int. 
  # Values ({0.0 ... 0.009} == 0) in rep() function of R.
  z <- rep(0.01, n_i)  
  
  for(i in 1:n_i) {
    x_new <- x_receptors[i]
    y_new <- y_receptors[i]
    z_new <- 0
    for(j in 1:n_j) {
      weight_j <- alpha_weights[j]
      
      x_j <- x_observations[j]
      y_j <- y_observations[j]
      
      r <- sqrt( ((x_j-x_new)^2) + ((y_j-y_new)^2) )
      m_distance <- mqrbf.calc_distance(r, smoothing_factor)
      
      z_new <- z_new + (weight_j * m_distance)
    }
    z[i] <- z_new
  }
  
  return(z)
  
}

mqrbf.calc_new_rsme <- function(weight, smoothing_factor = 78.234334){
  x_observations <- sf_observations$x
  y_observations <- sf_observations$y
  z_observations <- sf_observations$windspeed
  training_weight <- mqrbf.base_weight + weight
  
  n <- length(x_observations)
  residuals <- rep(0.01, n)
  
  for(i in 1:n) {
    x_train <- x_observations[-i]
    y_train <- y_observations[-i]
    z_train <- z_observations[-i]
    
    x_new <- as.vector(x_observations[i])
    y_new <- as.vector(y_observations[i])
    
    alpha_weights <- mqrbf.calc_alpha_weights(x_train, y_train, z_train, training_weight)
    
    z_new <- mqrbf.predict(
      x_new, y_new,
      x_train, y_train,
      alpha_weights, smoothing_factor
    )
    
    residuals[i] <- (z_observations[i] - z_new)
  }
  
  return(rmse(residuals))
  
}


## IDW
idw.base_weight <- 1.0

idw.calc_new_rsme <- function(weight, neighbors = 5){
  idp_weight <- idw.base_weight + weight
  
  model <- gstat(
    formula = windspeed ~ 1,
    data = sf_observations,
    nmax = neighbors,
    set = list(idp = idp_weight)
  )
  
  { sink(nullfile()); cv <- gstat.cv(model, nfold = nrow(sf_observations), verbose=FALSE); sink(); }
  
  return( rmse(cv$residual) )
}

trend_surface.best_fit_weights <- dynamic_model_train_vars(
  trend_surface.calc_new_rsme, 
  c("windspeed ~ x + y", "windspeed ~ x + y + I(x^2) + (x * y) + I(y^2) + I(x^3) + (I(x^2) * y) + (x * I(y^2)) + I(y^3)")
)

mqrbf.best_fit_weights <- dynamic_model_train_vars(
  mqrbf.calc_new_rsme,
  runif(10, min = min(sf_observations$x), max = max(sf_observations$x)),
  mqrbf.base_weight
)

idw.best_fit_weights <- dynamic_model_train_vars(
  idw.calc_new_rsme,
  seq(2, 8, 1),
  idw.base_weight
)



#cat(sprintf(
#"
#  trend_surface [linear] = (1 / windspeed) + %f
#  trend_surface [cubic] = (1 / windspeed) + %f
#  MQ-RBF = %f
#  IDW = %f
#", 
#  trend_surface.linear.best_fit_weight, trend_surface.cubic.best_fit_weight, 
#  mqrbf.best_fit_weight, idw.best_fit_weight
#))



