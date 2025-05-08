# This is to sample beta 
sample_beta <- function(y, current, y_hat, hyperparam) {
  n <- ncol(y)
  betas <- c()
  for (i in 1:n) {
    betas[i] <- sample_blr(y[, i], y_hat[, i] / current$beta[i],
                           mu = hyperparam$beta_prior_mu, sigma = hyperparam$beta_prior_sd,
                           V = diag(1, 1), a = 1, b = 1)$beta
  }
  return(betas)
}

##############################################################
# This is to sample blr
sample_blr <- function(y, X, mu, sigma, V, a, b) {
  n <- length(y)
  V_post <- solve(solve(V) + t(X) %*% X)
  mu_post <- V_post %*% (solve(V) %*% mu + t(X) %*% y)
  a_post <- a + n / 2
  b_post <- b + 1 / 2 * (mu %*% solve(V) %*% mu + t(y) %*% y - t(mu_post) %*% solve(V_post) %*% mu_post)
  
  beta <- mvtnorm::rmvt(n = 1, sigma = (b_post[1] / a_post) * V_post,
                        df = 2 * a_post, delta = mu_post)
  return(list(beta = beta,
              sigma2 = invgamma::rinvgamma(1, shape = a_post, rate = b_post)))
}



#########################################################
# This is to sample f 
sample_f <- function(y, theta, n_draws, nugget = 1e-6) {
  n_time <- nrow(y)
  n <- ncol(y)
  chain_f <- vector(mode = "list", length = n_draws)
  x <- seq(0, 1, length.out = n_time)
  K_f <- sq_exp_kernel(x, theta$rho, nugget = nugget)
  K_f_inv <- solve(K_f)
  Sigma_nu <- get_Sigma_nu(theta$phi, theta$sigma, n_time)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    A <- K_f_inv
    b <- matrix(0, n_time)
    # Sigma_i_inv <- diag(1 / theta$sigma^2, n_time)
    for (i in 1:n) {
      Sigma_y_i <- get_Sigma_y_i(theta$beta[i], K_f, Sigma_nu)
      K_i <- get_K_i22(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
      Sigma_i <- Sigma_y_i - t(K_i)%*%K_f_inv%*%K_i
      Sigma_i <- (Sigma_i + t(Sigma_i))/2
      
      Sigma_i_inv <- solve(Sigma_i)
      L <- K_i %*% K_f_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[, i] %*% G)
    }
    K_f_post <- solve(A)
    K_f_post <- (K_f_post + t(K_f_post))/2
    chain_f[[iter]] <- MASS::mvrnorm(n = 1, K_f_post %*% b, K_f_post)
  }
  f_draws <- matrix(unlist(lapply(chain_f, `[`)), nrow = n_time, ncol = n_draws)
  return(f_draws)
}


#####################################################################
# This is to sample f using the Nearest Neighbor Method
sample_f_NNGP <- function(y, theta, n_draws, m = 6, nugget = 1e-6) {
  n_time <- nrow(y)
  n <- ncol(y)
  chain_f <- vector(mode = "list", length = n_draws) # create a vector of 100 lists, each element is a list
  x <- seq(0, 1, length.out = n_time) # this part is to create a grid points between 0 and 1
  K_f <- sq_exp_kernel(x, theta$rho, nugget = nugget) 
  K_f_inv <- solve(K_f)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    A <- K_f_inv
    b <- matrix(0, n_time)
    Sigma_i_inv <- diag(1 / theta$sigma^2, m)
    for (i in 1:n) {
      K_i_result <- get_K_i_neighbor(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]),m=6)
      K_i <- K_i_result$vec
      nearest_neighbor_index <- K_i_result$target_index
      L <- K_i %*% K_f_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[nearest_neighbor_index, i] %*% G)
    }
    K_f_post <- solve(A)
    K_f_post <- (K_f_post + t(K_f_post))/2
    chain_f[[iter]] <- MASS::mvrnorm(n = 1, K_f_post %*% b, K_f_post)
  }
  f_draws <- matrix(unlist(lapply(chain_f, `[`)), nrow = n_time, ncol = n_draws)
  return(f_draws)
}



#########################################################################
# This is to sample f NNGP by sampling only one observation. 
sample_f_one_obs_NNGP <- function(y, theta, n_draws, m = 6, nugget = 1e-6) {
  n_time <- nrow(y)
  n <- 1
  chain_f <- vector(mode = "list", length = n_draws) # create a vector of 100 lists, each element is a list
  x <- seq(0, 1, length.out = n_time) # this part is to create a grid points between 0 and 1
  K_f <- sq_exp_kernel(x, theta$rho, nugget = nugget) 
  K_f_inv <- solve(K_f)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    A <- K_f_inv
    b <- matrix(0, n_time)
    Sigma_i_inv <- diag(1 / theta$sigma^2, m)
    for (i in 1:n) {
      K_i_result <- get_K_i_neighbor(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]),m=6)
      K_i <- K_i_result$vec
      nearest_neighbor_index <- K_i_result$target_index
      L <- K_i %*% K_f_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[nearest_neighbor_index, i] %*% G)
    }
    K_f_post <- solve(A)
    K_f_post <- (K_f_post + t(K_f_post))/2
    chain_f[[iter]] <- MASS::mvrnorm(n = 1, K_f_post %*% b, K_f_post)
  }
  f_draws <- matrix(unlist(lapply(chain_f, `[`)), nrow = n_time, ncol = n_draws)
  return(f_draws)
}


####################################################################
# This is to sample f with only one observation with the projection method
sample_f_one_obs_Projection <- function(y, theta, n_draws, m = 15,nugget = 1e-6) {
  n_time <- nrow(y)
  n <- 1
  chain_f <- vector(mode = "list", length = n_draws) # create a vector of 100 lists, each element is a list
  x <- seq(0, 1, length.out = n_time) # this part is to create a grid points between 0 and 1
  K_f <- sq_exp_kernel(x, theta$rho, nugget = nugget) 
  K_f_inv <- solve(K_f)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    A <- K_f_inv
    b <- matrix(0, n_time)
    Sigma_i_inv <- diag(1 / theta$sigma^2, m)
    for (i in 1:n) {
      K_i_result <- get_K_i_projection(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]),m=m)
      K_i <- K_i_result$vec
      zstar_index <- K_i_result$target_index
      L <- K_i %*% K_f_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[zstar_index, i] %*% G)
    }
    K_f_post <- solve(A)
    K_f_post <- (K_f_post + t(K_f_post))/2
    chain_f[[iter]] <- MASS::mvrnorm(n = 1, K_f_post %*% b, K_f_post)
  }
  f_draws <- matrix(unlist(lapply(chain_f, `[`)), nrow = n_time, ncol = n_draws)
  return(f_draws)
}



###########################################################################
# This is to sample f with projection method 
sample_f_projection <- function(y, theta, n_draws, m = 15,nugget = 1e-6) {
  n_time <- nrow(y)
  n <- ncol(y)
  chain_f <- vector(mode = "list", length = n_draws) # create a vector of 100 lists, each element is a list
  x <- seq(0, 1, length.out = n_time) # this part is to create a grid points between 0 and 1
  K_f <- sq_exp_kernel(x, theta$rho, nugget = nugget) 
  K_f_inv <- solve(K_f)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    A <- K_f_inv
    b <- matrix(0, n_time)
    Sigma_i_inv <- diag(1 / theta$sigma^2, m)
    for (i in 1:n) {
      K_i_result <- get_K_i_projection(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]),m=m)
      K_i <- K_i_result$vec
      zstar_index <- K_i_result$target_index
      L <- K_i %*% K_f_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[zstar_index, i] %*% G)
    }
    K_f_post <- solve(A)
    K_f_post <- (K_f_post + t(K_f_post))/2
    chain_f[[iter]] <- MASS::mvrnorm(n = 1, K_f_post %*% b, K_f_post)
  }
  f_draws <- matrix(unlist(lapply(chain_f, `[`)), nrow = n_time, ncol = n_draws)
  return(f_draws)
}

##########################################################################
# This is to sample f with updated projection method 
sample_f_projection2 <- function(y, theta, n_draws, m, nugget = 1e-6) {
  n_time <- nrow(y)
  n <- ncol(y)
  inducing_length <- m 
  chain_f <- vector(mode = "list", length = n_draws)
  chain_f_higher <- vector(mode="list",length=n_draws)
  x <- seq(0, 1, length.out = n_time)
  
  # choose the inducing points
  inducing_points_info <- get_inducing_points(x=x,m=m)
  inducing_points <- inducing_points_info$points
  inducing_indices <- inducing_points_info$indices
  
  # 
  K_z <- sq_exp_kernel2(inducing_points,inducing_points,theta$rho, nugget = nugget)
  K_z_inv <- solve(K_z)
  Sigma_nu <- get_Sigma_nu(theta$phi, theta$sigma, inducing_length)
  
  #
  K_x_z <- sq_exp_kernel2(x,inducing_points,theta$rho)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    
    A <- K_z_inv
    b <- matrix(0, inducing_length)
    
    for (i in 1:n) {
      Sigma_y_i <- get_Sigma_y_i(theta$beta[i], K_z, Sigma_nu)
      K_i <- get_K_i22(inducing_points, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
      Sigma_i <- Sigma_y_i - t(K_i)%*%K_z_inv%*%K_i
      Sigma_i <- (Sigma_i + t(Sigma_i))/2
      Sigma_i_inv <- solve(Sigma_i)
      
      L <- K_i %*% K_z_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[inducing_indices, i] %*% G) # this part is to loop over the entire Big Sigma Notation part
    }
    
    # This is K_f before the projection (still on the lower space)
    K_f_post_lower <- solve(A)
    K_f_post_lower <- (K_f_post_lower + t(K_f_post_lower))/2 # covariance before the projection 
    
    mu_post_lower <- K_f_post_lower %*% b    # mu before the projection 
    
    # Ensure cur_f is treated as a column vector
    cur_f <- as.matrix(MASS::mvrnorm(n = 1, mu_post_lower, K_f_post_lower))
    #cur_f <- t(cur_f)  # Ensure it's a column vector
    
    # Store the current sample in chain_f
    chain_f[[iter]] <- cur_f
    higher_f <- K_x_z %*% K_z_inv %*% cur_f
    chain_f_higher[[iter]] <- higher_f
  }
  
  f_draws <- matrix(unlist(lapply(chain_f_higher, `[`)), nrow = n_time, ncol = n_draws)
  
  return(f_draws)
}


########################################################################
# This is to sample f with updated Sigma_nu function 
sample_f_projection3 <- function(y, theta, n_draws, m, nugget = 1e-6) {
  n_time <- nrow(y)
  n <- ncol(y)
  inducing_length <- m 
  chain_f <- vector(mode = "list", length = n_draws)
  chain_f_higher <- vector(mode="list",length=n_draws)
  x <- seq(0, 1, length.out = n_time)
  
  # choose the inducing points
  inducing_points_info <- get_inducing_points(x=x,m=m)
  inducing_points <- inducing_points_info$points
  inducing_indices <- inducing_points_info$indices
  
  # 
  K_z <- sq_exp_kernel2(inducing_points,inducing_points,theta$rho, nugget = nugget)
  K_z_inv <- solve(K_z)
  Sigma_nu <- extract_sigma_nu(theta$phi, theta$sigma, n_time,inducing_length)
  
  #
  K_x_z <- sq_exp_kernel2(x,inducing_points,theta$rho)
  
  for (iter in 1:n_draws) {
    if (iter %% 10 == 0) cat(iter / 10)
    
    A <- K_z_inv
    b <- matrix(0, inducing_length)
    
    for (i in 1:n) {
      Sigma_y_i <- get_Sigma_y_i(theta$beta[i], K_z, Sigma_nu)
      K_i <- get_K_i22(inducing_points, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
      Sigma_i <- Sigma_y_i - t(K_i)%*%K_z_inv%*%K_i
      Sigma_i <- (Sigma_i + t(Sigma_i))/2
      Sigma_i_inv <- solve(Sigma_i)
      
      L <- K_i %*% K_z_inv
      G <- Sigma_i_inv %*% L
      A <- A + t(L) %*% G
      b <- b + t(y[inducing_indices, i] %*% G) # this part is to loop over the entire Big Sigma Notation part
    }
    
    # This is K_f before the projection (still on the lower space)
    K_f_post_lower <- solve(A)
    K_f_post_lower <- (K_f_post_lower + t(K_f_post_lower))/2 # covariance before the projection 
    
    mu_post_lower <- K_f_post_lower %*% b    # mu before the projection 
    
    # Ensure cur_f is treated as a column vector
    cur_f <- as.matrix(MASS::mvrnorm(n = 1, mu_post_lower, K_f_post_lower))
    #cur_f <- t(cur_f)  # Ensure it's a column vector
    
    # Store the current sample in chain_f
    chain_f[[iter]] <- cur_f
    higher_f <- K_x_z %*% K_z_inv %*% cur_f
    chain_f_higher[[iter]] <- higher_f
  }
  
  f_draws <- matrix(unlist(lapply(chain_f_higher, `[`)), nrow = n_time, ncol = n_draws)
  
  return(f_draws)
}



###############################################################################
# sample rho
sample_rho <- function(y, f, current, hyperparam) {
  
  n_time <- nrow(y)
  x <- seq(0, 1, length.out=n_time)
  
  proposed <- current
  proposed$rho <- propose_rho(current$rho, hyperparam$rho_proposal_sd)
  K_f_curr = sq_exp_kernel(x, current$rho, nugget = 1e-6)
  K_f_curr_inv = solve(K_f_curr)
  
  K_f_prop = sq_exp_kernel(x, proposed$rho, nugget = 1e-6)
  K_f_prop_inv = solve(K_f_prop)
  
  lik_current <- likelihood(y, f, current, K_f_curr, K_f_curr_inv); lik_current
  prior_current <- prior$rho(current$rho, hyperparam$rho_prior_shape, hyperparam$rho_prior_scale); 
  lik_proposed <- likelihood(y, f, proposed, K_f_prop, K_f_prop_inv); lik_proposed
  prior_proposed <- prior$rho(proposed$rho, hyperparam$rho_prior_shape, hyperparam$rho_prior_scale)
  prob <- exp(lik_proposed + prior_proposed - lik_current - prior_current)
  # cat("prob", prob, "\n")
  if (prob > runif(1)) {
    return(proposed$rho)
  } else {
    return(current$rho)
  }
}


###########################################################################
# sample tau
sample_tau <- function(y, f, current, hyperparam, K_f, K_f_inv) {
  proposed <- current
  proposed$tau <- propose_tau(current$tau, hyperparam$tau_proposal_sd)
  
  lik_current <- likelihood(y, f, current, K_f, K_f_inv)
  prior_current <- prior$tau(current$tau, hyperparam$tau_prior_sd)
  
  lik_proposed <- likelihood(y, f, proposed, K_f, K_f_inv)
  prior_proposed <- prior$tau(proposed$tau, hyperparam$tau_prior_sd)
  
  prob <- exp(lik_proposed + prior_proposed - lik_current - prior_current)
  if (prob > runif(1)) {
    return(proposed$tau)
  } else {
    return(current$tau)
  }
}


###########################################################################



