# This is the proposal function for rho 
propose_rho <- function(rho, rho_proposal_sd) {
  proposal <- rnorm(1, rho, rho_proposal_sd)
  return(proposal)
}


# This is the proposal function for tau
propose_tau <- function(tau, tau_proposal_sd) {
  proposal <- rep(NA, length(tau))
  n <- length(tau)
  Sigma <- tau_proposal_sd^2 * (diag(1, n))
  proposal[1:n] <- MASS::mvrnorm(n = 1, tau[1:n], Sigma)
  return(proposal)
}


# 

