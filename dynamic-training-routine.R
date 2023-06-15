# ====================== SIMULATION DATA LOAD (START) ======================

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
  sf_observations <- readRDS("Data/dynamic-training-routine.rds")
  
  # Split x and y coordinate from observations and receptors as separate columns
  sf_observations <- sf_observations %>%
    mutate(
      x = st_coordinates(.)[,1],
      y = st_coordinates(.)[,2]
    ) %>%
    select(x, y, windspeed)
  
  nrows_observations <- nrow(sf_observations)

# ====================== SIMULATION DATA LOAD (END) ======================

# ====================== SIMULATION SETTINGS (START) ======================
  
  settings.loss_function.default_power <- 0.68732
  
  settings.model_opt.polynomial_formulas <- c(
    "windspeed ~ x + y", 
    "windspeed ~ x + y + I(x^2) + (x * y) + I(y^2) + I(x^3) + (I(x^2) * y) + (x * I(y^2)) + I(y^3)",
    "windspeed ~ abs(x + y + (x^2) + (x * y) + (y^2) + (x^3) + ((x^2) * y) + (x * (y^2)) + (y^3))"
  )
  
  settings.model_opt.mqrbf_smoothing_factors <- runif(5, min=0.5, max=100.0)
  settings.model_opt.neighbors <- seq(5, 8, 1)
  
# ====================== SIMULATION SETTINGS (END) ======================
  
# ====================== SIMULATION HELPER FUNCTIONS (START) ======================
  
  # Define the function used to calculate training weights based on RMSE gain/loss
  loss_function <- function(rmse_old, rmse_new, power=settings.loss_function.default_power){
    rmse_diff <- abs(rmse_old - rmse_new)
    rmse_loss <- -rmse_diff^power
    return( abs(rmse_loss) )
  }
  
  # Define the function used to calculate the RMSE from residuals
  rmse <- function(residuals) {
    return(sqrt(sum((residuals)^2)/length(residuals)))
  }
  
  # Define weight functions
  weight_func.one_base_sub_halve <- function(
    old_weight, loss_value, weight_func_opts = FALSE
  ) {
    if( isFALSE(old_weight) ){ return( 1.0 + loss_value )  }
    
    if( loss_value < 0.0 ){
      weight_factor <- (abs(loss_value) / 2)
      if (old_weight > weight_factor){ return( old_weight - weight_factor ) }
      return (old_weight + weight_factor)
    } else{
      return( old_weight + loss_value ) 
    }
  }
  
  weight_func.no_base_add_exp_sub_halve <- function(
    old_weight, loss_value, exp = 1.67328
  ) {
    
    if( isFALSE(old_weight) ){ return( loss_value )  }
    
    if( loss_value < 0.0 ){
      weight_factor <- (abs(loss_value) / 2)
      if (old_weight > weight_factor){ return( old_weight - weight_factor ) }
      return (old_weight + weight_factor)
    } else{
      return( old_weight + (loss_value^exp ) ) 
    }
  }
  
  weight_func.rand_base_add_sub_inc <- function(
    old_weights, loss_value, weight_func_opts
  ){
    if(isFALSE(weight_func_opts)){ weight_func_opts <- list( n_weights = nrows_observations, alpha = 0.08232 ) }
    
    if( isFALSE(old_weights) ){ return(runif(
      weight_func_opts$n_weights, 0+loss_value+weight_func_opts$alpha, 1-loss_value-weight_func_opts$alpha) 
      )  
    }
    
    if( loss_value < 0.0 ){
      lowest_weight <- min(old_weights)
      if ((lowest_weight - loss_value) < 0.0){
        return( old_weights - lowest_weight + weight_func_opts$alpha )
      }
      return( old_weights - lowest_weight )
    } else{
      return( old_weights + loss_value ) 
    }
  }
  
  weight_func.variogram_base_add_exp <- function(
    old_vgm_weights, loss_value, weight_func_opts
  ) {
    if(isFALSE(weight_func_opts)){ weight_func_opts <- list( sill_power = 100, range_power = 10000, exp = 1.327347 ) }
    
    if( isFALSE(old_vgm_weights) ){
      
      new_weight <- weight_func.no_base_add_exp_sub_halve(FALSE, loss_value, weight_func_opts$exp)
      new_vgs_weights <- list(
        weight = new_weight, 
        psill = new_weight * weight_func_opts$sill_power, 
        range = new_weight * weight_func_opts$range_power
      )
      return(new_vgs_weights)
      
    }
    
    new_weight <- weight_func.no_base_add_exp_sub_halve(old_vgm_weights$weight, loss_value, weight_func_opts$exp)
    new_vgs_weights <- list(
      weight = new_weight, 
      psill = new_weight * weight_func_opts$sill_power, 
      range = new_weight * weight_func_opts$range_power
    )
    return(new_vgs_weights)
    
  }
  
# ====================== SIMULATION HELPER FUNCTIONS (END) ======================  
  
# ====================== SIMULATION FUNCTIONS (START) ======================
  
  # Define the function used to dynamically train a single variation of a model
  dynamic_model.train <- function(
    model_opt, model_opt_id,
    rmse_func, weight_func, weight_func_opts = FALSE,
    loss_function_power = settings.loss_function.default_power, 
    max_iterations = 20
  ){
    
    # Pre-calculate the first two RSME values.
    old_rsme <- rmse_func(weight_func(FALSE, 0.001, weight_func_opts), model_opt)
    new_rsme <- rmse_func(weight_func(FALSE, 0.005, weight_func_opts), model_opt) 
    
    # Pre-calculate the first loss_value
    loss_value <- loss_function(old_rsme, new_rsme, loss_function_power)
    
    # Save states of pre-calculated values
    weights <- weight_func(FALSE, loss_value, weight_func_opts) 
    
    training_results <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(training_results) <- c("model_opt", "train_iteration", "old_rsme", "new_rsme", "loss_value", "weights")
    
    old_rsme <- new_rsme
    
    for( i in 1:max_iterations ) {
      new_rsme <- rmse_func(weights, model_opt)
      loss_value <- loss_function(old_rsme, new_rsme, loss_function_power)
      
      if( new_rsme > old_rsme ) { 
        weights <- weight_func(weights, -loss_value, weight_func_opts)
      } else {
        weights <- weight_func(weights, loss_value, weight_func_opts) 
      }
      
      training_results[i,] <- c(model_opt_id, i, old_rsme, new_rsme, loss_value, toString(weights))
      old_rsme <- new_rsme
    }
    
    return(training_results)
    
  }
  
  # Define the function used to dynamically train all variations of a model
  dynamic_model.train_opts <- function(
    model_opts,
    rmse_func, weight_func, weight_func_opts = FALSE,
    loss_function_power = settings.loss_function.default_power,
    max_iterations = 20
  ){
    n_opts <- length(model_opts)
    if (n_opts < 1) { return(NULL) }
    
    # Pre-calculate first entry of opts
    model_opt <- model_opts[1]
    
    training_results <- dynamic_model.train(
      model_opt, 1, 
      rmse_func, weight_func, weight_func_opts,
      loss_function_power, max_iterations
    )
    
    if (n_opts == 1){ return(training_results) }
    
    for (i in 2:n_opts) {
      model_opt <- model_opts[i]
      
      training_result <- dynamic_model.train(
        model_opt, i, 
        rmse_func, weight_func, weight_func_opts,
        loss_function_power, max_iterations
      )
      
      training_results <- rbind(training_results, training_result)
    }
    
    return(training_results)
  }
  
# ====================== SIMULATION FUNCTIONS (END) ======================
  
# ====================== SIMULATION MODEL HELPER FUNCTIONS (START) ======================  
  
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
  
  od_kriging.vgm <- variogram(windspeed ~ 1, data = sf_observations)
  od_kriging.vgm.fit <- suppressWarnings(fit.variogram(od_kriging.vgm, vgm("Sph")))
  u_kriging.vgm <- variogram(windspeed ~ x + y, data = sf_observations)
  u_kriging.vgm.fit <- suppressWarnings(fit.variogram(u_kriging.vgm, vgm("Sph")))
  
# ====================== SIMULATION MODEL HELPER FUNCTIONS (END) ======================  
  
# ====================== SIMULATION MODEL FUNCTIONS (START) ======================

  rmse_func.trend_surface <- function(case_weight, formula){
    
    model <- train(
      as.formula(formula),
      method = "lm",
      data = sf_observations,
      trControl = trainControl(method = "LOOCV"),
      weights = (1 / windspeed) + case_weight
    )
    
    return( model$results$RMSE )
    
  }
  
  rmse_func.mqrbf <- function(alpha_seed, smoothing_factor){
    x_observations <- sf_observations$x
    y_observations <- sf_observations$y
    z_observations <- sf_observations$windspeed
    
    n <- length(x_observations)
    residuals <- rep(0.01, n)
    
    for(i in 1:n) {
      x_train <- x_observations[-i]
      y_train <- y_observations[-i]
      z_train <- z_observations[-i]
      
      x_new <- as.vector(x_observations[i])
      y_new <- as.vector(y_observations[i])
      
      alpha_weights <- mqrbf.calc_alpha_weights(x_train, y_train, z_train, alpha_seed)
      
      z_new <- mqrbf.predict(
        x_new, y_new,
        x_train, y_train,
        alpha_weights, smoothing_factor
      )
      
      residuals[i] <- (z_observations[i] - z_new)
    }
    
    return(rmse(residuals))
    
  }
  
  rmse_func.idw <- function(idp_power, neighbors){
    
    model <- gstat(
      formula = windspeed ~ 1,
      data = sf_observations,
      nmax = neighbors,
      set = list(idp = idp_power)
    )
    
    { sink(nullfile()); cv <- gstat.cv(model, nfold = nrow(sf_observations), verbose=FALSE); sink(); }
    
    return( rmse(cv$residual) )
  }
  
  rmse_func.od_kriging <- function(vgm_weights, neighbors){
    
    base_psill <-  od_kriging.vgm.fit$psill[1]
    base_range <- od_kriging.vgm.fit$range[2]
    new_sill <- base_psill + vgm_weights$psill
    new_range <- base_range + vgm_weights$range
    
    vgm_model <- suppressWarnings(fit.variogram(
      od_kriging.vgm, 
      vgm(
        psill = new_sill, 
        model = "Sph", 
        range = new_range
      )
    ))
    
    model <- suppressWarnings(gstat(
      formula = windspeed ~ 1,
      data = sf_observations,
      nmax = neighbors,
      model = vgm_model,
    ))
      
    { sink(nullfile()); cv <- gstat.cv(model, nfold = nrow(sf_observations), verbose=FALSE); sink(); }
    
    return( rmse(cv$residual) )
  }
  
  rmse_func.u_kriging <- function(vgm_weights, model_opts) {
    
    neighbors <- model_opts[[1]][1]
    formula <- model_opts[[1]][2]
    
    #print(sprintf("neighbors=%s formula=%s",neighbors,formula))
    
    base_psill <-  u_kriging.vgm.fit$psill[1]
    base_range <- u_kriging.vgm.fit$range[2]
    new_sill <- base_psill + vgm_weights$psill
    new_range <- base_range + vgm_weights$range
    
    vgm_model <- suppressWarnings(fit.variogram(
      u_kriging.vgm, 
      vgm(
        psill = new_sill, 
        model = "Sph", 
        range = new_range
      )
    ))
    
    model <- suppressWarnings(gstat(
      formula = as.formula(formula),
      data = sf_observations,
      nmax = neighbors,
      model = vgm_model,
    ))
    
    { sink(nullfile()); cv <- gstat.cv(model, nfold = nrow(sf_observations), verbose=FALSE); sink(); }
    
    return( rmse(cv$residual) )
  }

# ====================== SIMULATION MODEL FUNCTIONS (END) ======================

# ====================== SIMULATION DATA SAVE (START) ======================

  # Calculate training results for each model. 
  # training_results.trend_surface <- dynamic_model.train_opts(
  #   model_opts = settings.model_opt.polynomial_formulas,
  #   rmse_func = rmse_func.trend_surface,
  #   weight_func = weight_func.no_base_add_exp_sub_halve,
  #   max_iterations = 10
  # )
  # 
  # training_results.mqrbf <- dynamic_model.train_opts(
  #   model_opts = settings.model_opt.mqrbf_smoothing_factors,
  #   rmse_func = rmse_func.mqrbf,
  #   weight_func = weight_func.one_base_sub_halve,
  #   max_iterations = 10
  # )
  # 
  # training_results.idw <- dynamic_model.train_opts(
  #  model_opts = settings.model_opt.neighbors,
  #  rmse_func = rmse_func.idw,
  #  weight_func = weight_func.one_base_sub_halve,
  #  max_iterations = 10
  # )
  
  training_results.od_kriging <- dynamic_model.train_opts(
    model_opts = settings.model_opt.neighbors[c(1)],
    rmse_func = rmse_func.od_kriging,
    weight_func = weight_func.variogram_base_add_exp,
    weight_func_opts = list( sill_power = 100, range_power = 10000, exp=1.527347 ),
    max_iterations = 5
  )
  
  training_results.u_kriging <- dynamic_model.train_opts(
    model_opts = list(
      c(5, settings.model_opt.polynomial_formulas[1]),
      c(5, settings.model_opt.polynomial_formulas[3])
    ),
    rmse_func = rmse_func.u_kriging,
    weight_func = weight_func.variogram_base_add_exp,
    weight_func_opts = list( sill_power = 1000, range_power = 100000, exp=1.527347 ),
    max_iterations = 5
  )
  

# ====================== SIMULATION DATA SAVE (END) ======================





