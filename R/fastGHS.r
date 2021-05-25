#' Perform GHS with ECM/ICM
#' 
#' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the graphical horseshoe
#' 
#' @param X \eqn{n} by \eqn{p} matrix of data
#' @param theta initial value of the \eqn{p} by \eqn{p} precision matrix
#' @param sigma initial value of the \eqn{p} by \eqn{p} covariance matrix
#' @param Lambda_sq initial value of matrix of squared local shrinkage parameters
#' @param tau_sq initial value of squared global shrinkage parameter. If variables are grouped, a \eqn{ngroup} by \eqn{ngroup} matrix
#' @param method the method to use. Default is \eqn{ECM}. Other options include \eqn{ICM}
#' @param epsilon tolerance for the convergence assessment
#' @param maxitr maximum number of iterations
#' @param verbose logical indicator of printing information at each iteration
#' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for p<200
#' @param group grouping information, if variables are grouped. If provided, a vector of length \eqn{p} giving the group of each variable as a string, character or number. 
#' @param save_Q should the value of the objective function at each step be saved?
#' @return a fitted EMGS object
#' @export 
#' 
fastGHS <- function(X, theta=NULL,sigma=NULL,Lambda_sq=NULL, tau_sq = NULL, method= 'ECM',epsilon = 1e-5, maxitr = 1e5, verbose=TRUE, savepath = FALSE,  group=NULL, save_Q = F){

  p <- dim(X)[2]
  
  # Assign initial values, unless provided
  if(is.null(theta)){
    theta <- diag(1,p)
  }
  else{
    if(!isSymmetric(theta)){
      theta=as.matrix(Matrix::forceSymmetric(theta))
      cat('Initial theta not symmetric, forcing symmetry...')
    }
    if(ncol(theta)!= p | nrow(theta)!=p | !matrixcalc::is.positive.definite(theta+0)){
      cat('Error: initial theta must be pxp and positive definite \n')
      return()
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
    if(ncol(sigma)!= p | nrow(sigma)!=p | !matrixcalc::is.positive.definite(sigma+0)){
      cat('Error: initial sigma must be pxp and positive definite \n')
      return()
    } 
  }
  if(is.null(Lambda_sq)){
    Lambda_sq <- matrix(rep(1,p^2),ncol=p)
  }
  else{
    if(any(Lambda_sq<0)){
      cat('Error: negative Lambda_sq values are not allowed \n')
      return()
    } 
  }

  n <- dim(X)[1]
  S <- t(X) %*% X  #* n # Should we really multiply by n??

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
        cat('Error: negative tau_sq is not allowed \n')
        return()
      }
    }
    Tau_sq = S # Dummy variable
  }
  else{
    if(length(group)!=ncol(X)){
      cat('Error: number of group assignments and variables must match \n')
      return()
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
        cat('Error: negative tau_sq is not allowed \n')
        return()
      }
      Tau_sq = tau_sq;
    }
    tau_sq = 1; # Dummy variable
  }
  if(method=='ECM'){
    out <- ECM_GHS(as.matrix(X), S, theta , sigma, Lambda_sq, epsilon, verbose, maxitr, savepath, exist.group, group, N_groups, save_Q,tau_sq, Tau_sq, use_ICM = FALSE)
  }
  else if(method=='ICM'){
    out <- ECM_GHS(as.matrix(X), S, theta , sigma, Lambda_sq, epsilon, verbose, maxitr, savepath, exist.group, group, N_groups, save_Q,tau_sq, Tau_sq, use_ICM = TRUE)
  }
  else {
    out <- NULL
  }
  
  # Save outputs
  out$epsilon = epsilon
  out$maxitr = maxitr
  out$group = group
  
  class(out) = "fastGHS"
  return(out)
}