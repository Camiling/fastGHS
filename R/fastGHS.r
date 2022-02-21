#' Perform GHS with ECM
#' 
#' This function performs expectation-conditional-maximation for the graphical horseshoe
#' 
#' @param X \eqn{n} by \eqn{p} matrix of data.
#' @param theta initial value of the \eqn{p} by \eqn{p} precision matrix.
#' @param sigma initial value of the \eqn{p} by \eqn{p} covariance matrix.
#' @param Lambda_sq initial value of matrix of squared local shrinkage parameters.
#' @param tau_sq initial value of squared global shrinkage parameter. If variables are grouped, an \eqn{ngroup} by \eqn{ngroup} matrix.
#' @param method the method to use. Default is \code{ECM}. Other options include \code{ICM}.
#' @param AIC_selection logical. Should the global shrinkage parameters be selected with AIC?
#' @param AIC_eps if \code{AIC_selection=TRUE}, the convergence tolerance for the AIC convergence assessment.
#' @param tau_sq_min if \code{AIC_selection=TRUE}, the smallest value of \code{tau_sq} to consider.  
#' @param tau_sq_stepsize if \code{AIC_selection=TRUE}, the step-size to use in the grid for \code{tau_sq}. Optional argument.
#' @param tau_sq_max if \code{AIC_selection=TRUE}, the largest value of \code{tau_sq} to consider.  
#' @param epsilon tolerance for the convergence assessment.
#' @param maxitr maximum number of iterations.
#' @param verbose logical indicator of printing information at each iteration.
#' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for \eqn{p<200}.
#' @param group grouping information, if variables are grouped. If provided, a vector of length \eqn{p} giving the group of each variable as a string, character or number. 
#' @param save_Q should the value of the objective function at each step be saved?
#' @param fix_tau logical. Should \code{tau_sq} be fixed?
#' @param GHS_like logical. Should the GHS-like penalty be used?
#' @param stop_underflow should underflow be avoided by never allowing doubles to be smaller than the machine precision?
#' @param weights If provided, a vector of length \eqn{n} to weigh the obervations with. 
#' @param savepath_tau logical indicator of saving the value of \code{tau_sq} at each iteration in the ECM algorithm.
#' @param stop_overflow should measures be taken to avoid overflow? Should be set to \code{TRUE} for networks expected to be dense.
#' @return a fitted \code{fastGHS} object
#' @export 
#' 
fastGHS <- function(X, theta=NULL,sigma=NULL,Lambda_sq=NULL, tau_sq = NULL, method= 'ECM', AIC_selection=FALSE, AIC_eps = 1e-1, tau_sq_min =1e-4, tau_sq_stepsize= NULL, tau_sq_max = 20,
                    epsilon = 1e-5, maxitr = 1e4, verbose=TRUE, savepath = FALSE,  group=NULL, save_Q = F, fix_tau = FALSE, GHS_like = FALSE, stop_underflow = FALSE, weights=NA, savepath_tau=FALSE, stop_overflow=FALSE){

  # If the GHS-like penalty is used, tau_sq represents the fixed shrinkage parameter 'a' and N is the matrix of nu_ij's
  
  if(GHS_like){
    fix_tau = TRUE # a is not to be updated
    group = NULL # Groups are not implemented for GHS-like penalty
    use_ICM = FALSE
    save_Q = FALSE
    # random starting point
    theta = matrix(runif(ncol(X)^2,1e-3,0.3),nrow=ncol(X))
    diag(theta) = 1
    theta = as.matrix(Matrix::nearPD(theta)$mat)
    sigma = as.matrix(Matrix::nearPD(solve(theta))$mat)
  }
  
  p <- dim(X)[2]
  
  # Assign initial values, unless provided
  if(is.null(theta)){
    theta <- diag(1,p)
  }
  else{
    theta=as.matrix(theta)
    if(!isSymmetric(theta)){
      theta=as.matrix(Matrix::forceSymmetric(theta))
      cat('Initial theta not symmetric, forcing symmetry...')
    }
    if(ncol(theta)!= p | nrow(theta)!=p | !matrixcalc::is.positive.definite(theta+0)){
      stop('Initial theta must be pxp and positive definite \n')
    } 
  }
  if(is.null(sigma)){
    sigma <- diag(1,p)
  }
  else{
    if(!isSymmetric(sigma)){
      sigma=as.matrix(Matrix::forceSymmetric(sigma))
      cat('Initial sigma not symmetric, forcing symmetry...\n')
    }
    if(ncol(sigma)!= p | nrow(sigma)!=p | !matrixcalc::is.positive.definite(as.matrix(sigma+0))){
      stop('Initial sigma must be pxp and positive definite \n')
    } 
  }
  if(is.null(Lambda_sq)){
    Lambda_sq <- matrix(rep(1,p^2),ncol=p)
  }
  else{
    if(any(Lambda_sq<0)){
      stop('Negative Lambda_sq values are not allowed \n')
    } 
  }

  n <- dim(X)[1]
  if(all(!is.na(weights))){
    if(round(sum(weights),5)!=1){
      stop('Weights must sum to one \n')
    }
    if(any(weights<0)){
      stop('Negative weights not allowed \n')
    }
    if(length(weights)!=n){
      stop('Number of weights must match the number of samples \n')
    }
    X.w = X*sqrt(weights)
    S <- (n-1)/(1-sum(weights^2)) * (t(X.w) %*% X.w)
  }
  else{
     S <- t(X) %*% X
  }
  if(GHS_like){
    S <- S/n
  }

  # Check if vars should be grouped
  if(is.null(group)){
    group <- rep(0, dim(X)[2])
    exist.group <- 0
    N_groups <- matrix(rep(0,dim(X)[2]*2),nrow=2)
    if(is.null(tau_sq)){
      tau_sq <- 1
    }
    else{
      if(tau_sq<0){
        stop('Negative tau_sq is not allowed \n')
      }
    }
    Tau_sq = S # Dummy variable
  }
  else{
    if(length(group)!=ncol(X)){
      stop('Number of group assignments and variables must match \n')
    }
    group <- match(group, unique(group)) - 1
    exist.group <- length(unique(group))
    # Create matrix for storing number of observations in each group combination. 
    n_groups <- as.vector(table(group)) # Number of observations in each group
    N_groups_temp1 <- matrix(rep(n_groups, exist.group),ncol=exist.group,byrow=F)
    N_groups_temp2 <- matrix(rep(n_groups, exist.group),ncol=exist.group,byrow=T)
    diag(N_groups_temp2) <- 0
    N_groups <- N_groups_temp1 + N_groups_temp2
    
    # Create ngroup x ngroup matrix of between- and within-group shrinkage parameters. 
    if(is.null(tau_sq)){
      Tau_sq <- matrix(rep(1,exist.group^2),ncol=exist.group)
    }
    else{
      if(tau_sq<0){
        stop('Negative tau_sq is not allowed \n')
      }
      Tau_sq = tau_sq;
    }
    tau_sq = 1; # Dummy variable
  }
  machine_eps = .Machine$double.eps
  
  if(!method %in% c('ECM', 'ICM') ){
    stop('Method must be either ECM or ICM')
  }
  use_ICM = method=='ICM'
  
  if(AIC_selection){
    if(is.null(tau_sq_min)){
      tau_sq_min = 1e-3
    }
    if(is.null(tau_sq_stepsize)){
      tau_sq_vec = seq(tau_sq_min,tau_sq_max,by=2e-1)
    }
    else{
      tau_sq_vec = seq(tau_sq_min,tau_sq_max,by=tau_sq_stepsize)
    }
    out <- ECM_GHS_AIC(as.matrix(X), as.matrix(S), theta , sigma, Lambda_sq, AIC_eps, tau_sq_vec, epsilon, verbose, maxitr, savepath, exist.group, group, N_groups, save_Q,tau_sq, Tau_sq, machine_eps, use_ICM = use_ICM, GHS_like = GHS_like, stop_underflow=stop_underflow, stop_overflow=stop_overflow)
    out$AIC_scores = out$AIC_scores[1:out$ind]
    out$tau_sq_vals = tau_sq_vec[1:out$ind] 
  }
  else{
    out <- ECM_GHS(as.matrix(X), as.matrix(S), theta , sigma, Lambda_sq, epsilon, verbose, maxitr, savepath, exist.group, group, N_groups, save_Q,tau_sq, Tau_sq, machine_eps, use_ICM = use_ICM, fix_tau = fix_tau, GHS_like = GHS_like, stop_underflow=stop_underflow, savepath_tau=savepath_tau)
  }
  if(savepath_tau){
    out$tau_sq_path = out$tau_sq_path[1:out$count]
  }
  
  # Save outputs
  out$epsilon = epsilon
  out$maxitr = maxitr
  out$group = group
  class(out) = "fastGHS"
  return(out)
}