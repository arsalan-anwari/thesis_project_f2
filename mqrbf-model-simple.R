mqrbf <- function(data, alpha_seed = 1.0, smoothing_factor = 0.0) {
  n <- nrow(data)
  K <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- sqrt(
        1 + ( alpha_seed * 
                ( (sqrt( 
                  ( ( data$x[i] - data$x[j] )^2 ) + 
                    ( ( data$y[i] - data$y[j] )^2 ) 
                )^2) 
                ) 
        )
      )
      
    }
  }
  
  alpha_weights <- solve(K) %*% data$windspeed
  alpha_weights <- as.vector(alpha_weights)
  
  return(list(
    data = data, 
    alpha_weights = alpha_weights, 
    alpha_seed = alpha_seed, 
    smoothing_factor = smoothing_factor
  ))
}


mqrbf.distance_func <- function(r, c){
  return ( sqrt( (r^2) + (c^2) ) )
}

predict.mqrbf <- function(object, newdata){
  
  x_receptors <- newdata$x
  y_receptors <- newdata$y
  x_observations <- object$data$x
  y_observations <- object$data$y
  alpha_weights <- object$alpha_weights
  smoothing_factor <- object$smoothing_factor
  
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
      m_distance <- mqrbf.distance_func(r, smoothing_factor)
      
      z_new <- z_new + (weight_j * m_distance)
    }
    z[i] <- z_new
  }
  
  return(z)
  
}

train.mqrbf <- function(model){
  
  x_observations <- model$data$x
  y_observations <- model$data$y
  z_observations <- model$data$windspeed
  alpha_seed <- model$alpha_seed
  smoothing_factor <- model$smoothing_factor
  
  n <- length(x_observations)
  residuals <- rep(0.01, n)
  
  for(i in 1:n) {
    x_train <- x_observations[-i]
    y_train <- y_observations[-i]
    z_train <- z_observations[-i]
    data_train <- data.frame(x = x_train, y = y_train, windspeed = z_train)
    
    x_new <- as.vector(x_observations[i])
    y_new <- as.vector(y_observations[i])
    newdata <- list(x = x_new, y = y_new)
    
    model <- mqrbf(data_train, alpha_seed, smoothing_factor)
    
    z_new <- predict.mqrbf(model, newdata)
    
    residuals[i] <- (z_observations[i] - z_new)
  }
  
  return(list(residuals=residuals))
}