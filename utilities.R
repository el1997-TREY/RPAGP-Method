# This is the function to generate simulated data 
generate_data <- function(n, n_time, theta) {
  x <- seq(0, 1, length.out = n_time)
  K_f <- sq_exp_kernel(x, theta$rho, nugget = 1e-6)
  K_f_inv <- solve(K_f)
  
  f <- MASS::mvrnorm(n = 1, rep(0, n_time), K_f)
  y <- matrix(NA, nrow = n_time, ncol = n)
  z <- matrix(NA, nrow = n_time, ncol = n)
  mu <- matrix(NA, nrow = n_time, ncol = n)
  for (i in 1:n) {
    K_i <- get_K_i(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
    mu[, i] <- K_i %*% K_f_inv %*% f
    z[, i] <- arima.sim(model = list(ar = theta$phi), sd = theta$sigma, n = n_time)
    y[, i] <- mu[, i] + z[, i]
  }
  return(list(y = y, f = f, z = z, mu = mu))
}

###################################################################################
# This is the function to directly get the inducing points for the projection method.  
get_inducing_points <- function(x, m){
  # Set m (the number of intervals) 
  # Generate the sequence from 0 to 1
  seq_values <- x
  # Define the number of points in the sequence
  n <- length(seq_values)
  # Calculate the breakpoints for dividing the sequence into m intervals
  breakpoints <- seq(1, n, length.out = m + 1)
  # Initialize vectors to store sampled points and indices
  sampled_points <- numeric(m)
  sampled_indices <- numeric(m)
  # Loop through each interval, sample one point, and store its value and index
  for (i in 1:m) {
    # Get the start and end indices of the current interval
    start_idx <- ceiling(breakpoints[i])
    end_idx <- floor(breakpoints[i + 1])
    # Randomly sample one index from the interval
    sampled_indices[i] <- sample(start_idx:end_idx, 1)
    # Store the corresponding value
    sampled_points[i] <- seq_values[sampled_indices[i]]
  }
  return (list(points = sampled_points, 
               indices = sampled_indices))
}


###################################################################################
# This is the function to directly get the inducing points for the projection method.  
get_pseudo_inputs <- function(A, m_req,n_time) {
  # Inputs:
  # A: Initial RP partition scheme (list of layers with blocks)
  # X: Observed dataset (vector of points)
  # m_req: Minimum number of pseudo-inputs required per block
  
  # Initialize lists to store selected points and their indices
  selected_points <- list()
  selected_indices <- list()
  
  # Iterate through layers from bottom to top
  for (l in seq(length(A), 1, by = -1)) {
    if (length(A[[l]]) == 0) next  # Skip empty layers
    
    # Iterate through each block in the l-th layer
    for (j in seq_along(A[[l]])) {
      block <- A[[l]][[j]]
      
      # Find points and their indices in the current block
      block_indices <- which(X >= block$start & X <= block$end)
      block_points <- X[block_indices]
      
      # Check if the block meets the pseudo-input requirement
      if (length(block_points) >= m_req) {
        # Select m_req pseudo-input points
        if (length(block_points) >= m_req) {
          selected_indices_in_block <- sample(block_indices, m_req)
        } else {
          # If fewer points, select all available points
          selected_indices_in_block <- block_indices
        }
        
        # Store selected points and their indices
        selected_points[[length(selected_points) + 1]] <- X[selected_indices_in_block]
        selected_indices[[length(selected_indices) + 1]] <- selected_indices_in_block
      }
    }
  }
  
  # Return selected points and their indices
  return(list(points = unlist(selected_points), indices = unlist(selected_indices)))
}


################################################################################
# This is the function to get Cov(y_i, f)
get_K_i <- function(x, theta) {
  n_time <- length(x)
  K <- matrix(NA, n_time, n_time)
  for (i in 1:n_time) {
    for (j in 1:n_time) {
      K[i, j] <- exp(-theta$rho^2 / 2 * (x[i] - x[j] - theta$tau)^2)
    }
  }
  return(theta$beta * K)
}

######################################################################
# This is the function to get Cov(y_i,f) in the context of nearest neighbors 
get_K_i_neighbor <- function (x,theta,m){
  n_time <- length(x)
  K <- matrix(NA, m, n_time)
  for (i in 1:n_time){
    cur_time_point <- x[i]
    distances <- abs(x - cur_time_point)
    nearest_indices <- order(distances)[1:m]
    nearest_points <- x[nearest_indices]
    for (j in 1:m){
      K[j,i] <- exp(-theta$rho^2 / 2 * (x[i] - nearest_points[j] - theta$tau)^2)
    }
  }
  value_result <- theta$beta * K
  final_result <- list(vec = value_result, 
                       target_index = nearest_indices)
  return(final_result)
}

###################################################################
# This is the function to get Cov(y_i, f) in the context of projection methods 
get_K_i_projection <- function(x,theta,m){
  n_time <- length(x)
  K <- matrix(NA, m, n_time)
  
  z_star_idxs <- numeric(m) # get m 0's. 
  z_stars<- numeric(m)
  interval_size <- ceiling(n_time / m)
  for (i in 1:m) {
    # Define the start and end index of the current interval
    start_idx <- (i - 1) * interval_size + 1
    end_idx <- min(i * interval_size, n)  # Make sure not to exceed the length of the sequence
    # Sample one point from the current interval
    z_stars[i] <- sample(x[start_idx:end_idx], 1)
    z_star_idxs[i] <- which(x == z_stars[i])
  }
  
  for (i in 1:n_time){
    for (j in 1:m){
      K[j,i] <- exp(-theta$rho^2 / 2 * (x[i] - z_stars[j] - theta$tau[j])^2)
    }
  }
  value_result <- theta$beta * K
  final_result <- list(vec = value_result, 
                       target_index = z_star_idxs)
  return(final_result)
}


###################################################################
# This is the second way to get Cov(y_i,f)
get_Ki2 <- function(x,z,theta) {
  n_time <- length(x)
  reduced_time <- length(z)
  K <- matrix(NA, n_time, reduced_time)
  for (i in 1:n_time) {
    for (j in 1:reduced_time) {
      K[i, j] <- exp(-theta$rho^2 / 2 * (x[i] - x[j] - theta$tau)^2)
    }
  }
  return(theta$beta * K)
}


#################################################################
# This is the most efficient second way to get Cov(y_i,f)
get_K_i22 <- function(x, theta) {
  n_time <- length(x)
  
  # Compute the pairwise differences and adjust with theta$tau
  dist_sq <- (outer(x, x, "-") - theta$tau)^2  # Vectorized computation of (x[i] - x[j] - theta$tau)^2
  
  # Compute the kernel matrix
  K <- exp(-theta$rho^2 / 2 * dist_sq)
  
  # Multiply the entire matrix by theta$beta
  return(theta$beta * K)
}

##################################################################
# This is the function to get Sigma_nu
get_Sigma_nu <- function(phi, sigma, n_time)
{
  Sigma_nu_U <- matrix(0, n_time, n_time)
  tmp <- tacvfARMA(phi = phi, maxLag = n_time, sigma2 = sigma^2)
  for (tt in 1:n_time) {
    Sigma_nu_U[tt, tt:(n_time)] = tmp[1:(n_time-tt+1)]
  }
  Sigma_nu <- ultosymmetric(Sigma_nu_U)
  return(Sigma_nu)
}


##################################################################
# This is the first function to get Sigma_y_i
get_Sigma_y_i <- function(beta_i, K_f, Sigma_nu)
{
  return((beta_i^2)*K_f + Sigma_nu)
}


###################################################################
# This is the function to get y_hat 
get_y_hat <- function(i, f, theta, K_f_inv) {
  x <- seq(0, 1, length.out = length(f))
  K_i <- get_K_i22(x, list(tau = theta$tau[i], beta = theta$beta[i], rho = theta$rho))
  mu <- K_i %*% K_f_inv %*% f
  return(mu)
}


###################################################################
# This is the function to get y_hat matrix 
get_y_hat_matrix <- function(y, f, theta, K_f_inv) {
  y_hat <- matrix(nrow = nrow(y), ncol = ncol(y))
  
  for (i in 1:ncol(y_hat)) {
    y_hat[, i] <- get_y_hat(i, f, theta, K_f_inv)
  }
  return(y_hat)
}

######################################################################
# This is the function to get LLPPest
getLLPest <- function(singleTrialEstimates, categories)
{
  n_time <- dim(singleTrialEstimates)[1]
  n <- dim(singleTrialEstimates)[2]
  n_final <- dim(singleTrialEstimates)[3]
  T_LPP <- length(300:700)
  
  categories_labels <- c("High", "Low", "Neutral")
  out <- list()
  
  for (jj in 1:3) {
    idxs = which(categories == categories_labels[jj])
    n_jj <- length(idxs)
    out[[jj]] <- sapply(1:n_final, function(t) sum(singleTrialEstimates[, idxs, t]))/(T_LPP*n_jj)
  }
  
  df_out <- data.frame(iter = 1:n_final, High = out[[1]], Low = out[[2]], Neutral = out[[3]])
  df_out <- df_out %>% mutate(H_N = High - Neutral, 
                              L_N = Low - Neutral, 
                              H_L = High - Low)
}



#########################################################################
# This is the function to get modified LLPPtest 
getLPPest_mod <- function(singleTrialEstimates, categories)
{
  n_time <- dim(singleTrialEstimates)[1]
  n <- dim(singleTrialEstimates)[2]
  n_final <- dim(singleTrialEstimates)[3]
  T_LPP <- length(300:700)
  tmp <- apply(singleTrialEstimates, c(2, 3), sum)
  categories_labels <- c("High", "Low", "Neutral")
  
  
  out <- list()
  for (jj in 1:3) {
    idxs = which(categories == categories_labels[jj])
    n_jj <- length(idxs)
    out[[jj]] <- (apply(tmp[idxs, ], 2, sum)/(T_LPP*n_jj))
  }
  
  df_out <- data.frame(iter = 1:n_final, High = out[[1]], 
                       Low = out[[2]], Neutral = out[[3]])
  df_out <- df_out %>% mutate(H_N = High - Neutral, 
                              L_N = Low - Neutral, 
                              H_L = High - Low)
  return(df_out)
}


##########################################################################
# This is the first picture get 
getLPPest_orderedTrials <- function(singleTrialEstimates, categories, idxSortedTrials)
{
  n_time <- dim(singleTrialEstimates)[1]
  n <- dim(singleTrialEstimates)[2]
  n_final <- dim(singleTrialEstimates)[3]
  T_LPP <- length(300:700)
  categories_labels <- c("High", "Low", "Neutral")
  
  trial_categories <- rep(NA, n)
  medians <- numeric(n)
  lowers <- numeric(n)
  uppers <- numeric(n)
  loess <- numeric(n)
  loess_lwr <- numeric(n)
  loess_upr <- numeric(n)
  
  for (jj in 1:3) {
    idxs = which(categories == categories_labels[jj])
    idxs_aux = which(categories[idxSortedTrials] == categories_labels[jj])
    trial_categories[idxs_aux] = categories_labels[jj]
    
    n_jj <- length(idxs)
    tmp <- matrix(NA, n_jj, n_final)
    for (tt in 1:n_final) {
      tmp[, tt] <- apply(singleTrialEstimates[, idxs, tt], 2, sum)/T_LPP
    }
    quantiles_jj <- sapply(1:n_jj, function(i) {quantile(tmp[i, ], probs = c(0.025, 0.5, 0.975))})
    medians[idxs_aux] <- quantiles_jj[2, ]
    lowers[idxs_aux] <- quantiles_jj[1, ]
    uppers[idxs_aux] <- quantiles_jj[3, ]
    
    # l <- loess.sd(medians[idxs_aux] ~ idxs_aux, nsigma=1.96)
    
    lci <- loess.ci(medians[idxs_aux], idxs_aux, plot=FALSE, span=0.60)
    
    loess[idxs_aux] <- lci$loess
    loess_lwr[idxs_aux] <- lci$lci
    loess_upr[idxs_aux] <- lci$uci
  }
  
  df_out <- data.frame(trial = 1:n, 
                       category = trial_categories, 
                       med = medians, 
                       lwr = lowers, 
                       upr = uppers, 
                       loess = loess,
                       loess_lwr = loess_lwr, 
                       loess_upr = loess_upr)
  return(df_out)
}

#######################################################################
getLPPest_singleTrials <- function(singleTrialEstimates, categories)
{
  T_LPP <- length(300:700)
  tmp <- apply(singleTrialEstimates, c(2, 3), sum)/T_LPP
  df <- reshape2::melt(tmp, varnames = c("trial", "iter")) %>% 
    mutate(category = categories[trial])
  return(df)
}


######################################################################
getSigma_y_i_f <- function(i, x, theta, K_f, K_f_inv, Sigma_nu)
{
  Sigma_y_i <- get_Sigma_y_i(theta$beta[i], K_f, Sigma_nu)
  K_i <- get_K_i(x, list(rho = theta$rho, tau = theta$tau[i], beta = theta$beta[i]))
  Sigma_i <- Sigma_y_i - t(K_i)%*%K_f_inv%*%K_i
  Sigma_i <- (Sigma_i + t(Sigma_i))/2
  return(Sigma_i)
}


#######################################################################
getSingleTrialEstimates <- function(results, burn_in)
{
  n_iter <- length(results$chain)
  n_final <- length((burn_in+1):n_iter)
  
  # - beta 
  chain_beta_burned = matrix(NA, n, n_final)
  ss <- 1
  for (tt in (burn_in+1):n_iter) {
    chain_beta_burned[, ss] <- results$chain[[tt]]$beta 
    ss <- ss + 1
  }
  # - f 
  chain_f_burned = matrix(NA, n_time, n_final)
  ss <- 1
  for (tt in (burn_in+1):n_iter) {
    chain_f_burned[, ss] <- results$chain_f[[tt]]
    ss <- ss + 1
  }
  # - 
  y_hat = array(NA, dim = c(n_time, n, n_final))
  for (tt in 1:n_final) {
    for (ii in 1:n) {
      y_hat[, ii, tt] = chain_beta_burned[ii, tt] * chain_f_burned[, tt]
    }
  }
  return(y_hat)
}


###################################################################
getSortedTrials <- function(subj_id)
{
  dat_raw_ord <- read_delim(glue::glue("Only_WEIGHT_PicByPic_ERPsnocigPresentatioOrder.txt"), 
                            delim = "\t", col_names = T) 
  # filter categories of interest
  dat_ord <- tidyr::pivot_longer(dat_raw_ord, cols = 8:232) %>%
    filter(SubjectID == subj_id, CATOKnumFD %in% c("PH", "pl", "ul", "UH", "np")) %>% 
    rename(trial = prognum) %>% 
    mutate(category = recode_category(CATOKnumFD), 
           time = as.numeric(name)) %>% 
    dplyr::select(SubjectID, trial, category, time, value, PresentationOrder) %>% 
    filter(time >= time_lwr, time < time_upr) %>% 
    arrange(trial, PresentationOrder, time) 
  
  selected_trials <- unique(dat_ord %$% trial)
  n <- length(selected_trials)
  presentation_order <- numeric(n)
  
  for (i in 1:n) {
    presentation_order[i] <- dat_ord %>% filter(trial == selected_trials[i]) %$% PresentationOrder[1]
  }
  return(order(presentation_order))
}

#######################################################################
getSummaryOutput <- function(results, dat_trials, y, burn_in)
{
  n <- dim(y)[2]
  n_time <- dim(y)[1]
  n_iter <- length(results$chain)
  n_final <- length((burn_in+1):n_iter)
  
  y_hat <- getSingleTrialEstimates(results, burn_in)
  probs = c(0.025, 0.5, 0.975)
  
  y_hat_quantiles <- array(NA, dim= c(3, n_time, n))
  for (ii in 1:n) {
    y_hat_quantiles[, , ii] = sapply(1:n_time, function(t) quantile(y_hat[t, ii, ], 
                                                                    probs = probs))
  }
  lower <- reshape2::melt(y_hat_quantiles[1, , ], varnames = c("time", "trial"))$value
  median <- reshape2::melt(y_hat_quantiles[2, , ], varnames = c("time", "trial"))$value
  upper <- reshape2::melt(y_hat_quantiles[3, , ], varnames = c("time", "trial"))$value
  
  out <- dat_trials %>% mutate(lwr = lower, med = median, upr = upper)
  
  return(tibble(out))
}


####################################################################
NormalizeToRange <- function(x, r_min, r_max, t_min, t_max)
{
  out = ((x - r_min)/(r_max - r_min))*(t_max - t_min) + t_min
  return(out)
}


###################################################################
recode_category <- function(CATOKnumFD) {
  n <- length(CATOKnumFD)
  recoded <- rep(NA, n)
  for (i in 1:n) {
    if (CATOKnumFD[i] == "UH" | CATOKnumFD[i] == "PH") {
      recoded[i] <- "High"
    } else if (CATOKnumFD[i] == "ul" | CATOKnumFD[i] == "pl") {
      recoded[i] <- "Low"
    } else {
      recoded[i] <- "Neutral"
    }
  }
  return(recoded)
}

###################################################################
sample_AR <- function(z, ar.order)
{
  while(TRUE) {
    tmp <- suppressWarnings(gibbs_ar(c(z), 
                                     ar.order = ar.order, 
                                     Ntotal = 2, burnin =1))
    if (is.stationary(tmp$rho)) return(tmp)
  }
}

#####################################################################
sq_exp_kernel <- function(x, rho, alpha = 1, nugget = 0) {
  K <- toeplitz(alpha^2 * exp(-rho^2 / 2 * x^2))
  diag(K) <- diag(K) + nugget
  return(K)
}

###########################################################
sq_exp_kernel_projection <- function(x, rho, alpha = 1, nugget = 0) {
  n <- length(x)
  K <- matrix(0, n, n)
  # Compute the pairwise squared distances
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- alpha^2 * exp(-((x[i] - x[j])^2) / (2 * rho^2))
    }
  }
  # Add nugget to the diagonal to ensure numerical stability
  diag(K) <- diag(K) + nugget
  return(K)
}

###############################################################
sq_exp_kernel2 <- function(x1, x2, rho, alpha = 1, nugget = 0.0001) {
  # Precompute constants
  rho_sq_half <- rho^2 / 2
  # Vectorized computation of squared distances between x1 and x2
  dist_sq <- outer(x1, x2, "-")^2
  # Compute the kernel matrix using vectorized operations
  K <- alpha^2 * exp(-rho_sq_half * dist_sq)
  
  # Add the nugget to the diagonal if x1 and x2 are the same
  if (all(x1 == x2)) {
    diag(K) <- diag(K) + nugget
  }
  return(K)
}

#################################################################
sq_exp_kernel2_sparse <- function(x1, x2, rho, alpha = 1, thresh = 1e-5,nugget = 0.0001) {
  # Precompute constants
  rho_sq_half <- rho^2 / 2
  # Vectorized computation of squared distances between x1 and x2
  dist_sq <- outer(x1, x2, "-")^2
  # Compute the kernel matrix using vectorized operations
  K <- alpha^2 * exp(-rho_sq_half * dist_sq)
  
  K[K< thresh] = 0 
  # Add the nugget to the diagonal if x1 and x2 are the same
  if (all(x1 == x2)) {
    diag(K) <- diag(K) + nugget
  }
  return(K)
}

#################################################################
sq_exp_kernel21 <- function(x1,x2,rho,alpha=1,nugget=0.0001){
  K <- matrix(NA,nrow=length(x1),ncol=length(x2))
  for(i in 1:length(x1)){
    for (j in 1:length(x2)){
      K[i,j] <- alpha^2 * exp(-rho^2 / 2 * (x1[i] - x2[j])^2)
    }
    
    #Add the nugget to the diagonal if x1 and x2 are the same (x1 == x2)
    if (all(x1 == x2)) {
      diag(K) <- diag(K) + nugget
    }
  }
  return (K)
}


#################################################################
ultosymmetric <- function(m) 
{
  m = m + t(m) - diag(diag(m))
  return (m)
}


###########################################################################
extract_sigma_nu <- function(phi,sigma,n_time,m){
  Ori_Sigma_nu <- get_Sigma_nu(phi,sigma,n_time)
  indices <- seq(1,n_time,length.out=m)
  indices <- round(indices)
  target_matrix <- Ori_Sigma_nu[indices, indices]
  return(target_matrix)
}

########################################################################
# This function is used to calculate the MSE along with the two plots

Calculate_MSE <- function(true_f, melt_result){
  library(dplyr)
  median_estimate <- melt_result %>% 
    group_by(time) %>% 
    summarise(median = median(value)) %>%
    dplyr::select(median)
  MSE <- mean((true_f-unlist(median_estimate))^2)
  return(MSE)
}






