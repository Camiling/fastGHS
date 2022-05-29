//#include "ECM_GHS.h"
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

#include "E_xi.h"
#include "M_theta.h"
#include "M_tau.h"
#include "E_Nu.h"
#include "M_lambda.h"
#include "Q_val.h"

using namespace std;
using namespace Rcpp;
using namespace arma;



//' Perform GHS with ECM
//' 
//' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the graphical horseshoe
//' 
//' @param X \eqn{n} by \eqn{p} matrix of data
//' @param S \eqn{p} by \eqn{p} scatter matrix of the data
//' @param theta initial value of the \eqn{p} by \eqn{p} precision matrix
//' @param sigma initial value of the \eqn{p} by \eqn{p} covariance matrix
//' @param Lambda_sq initial value of matrix of squared local shrinkage parameters
//' @param epsilon tolerance for the convergence assessment
//' @param verbose logical indicator of printing information at each iteration
//' @param maxitr maximum number of iterations
//' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for p<200
//' @param save_Q should the value of the objective function at each step be saved?
//' @param tau_sq initial value of squared global shrinkage parameter. If exist_group==T, a dummy value should be provided
//' @param machine_eps numerical. The machine precision
//' @param fix_tau should tau_sq be fixed?
//' @param savepath_tau if fix_tau==F, should its path be saved?
//' @param stop_underflow Should tricks to avoid underflow be used?
//' 
//' @return A List with resulting ECM estimates, and saved path and objective function convergence information if requested
//' 
//' @export
// [[Rcpp::export]]
List ECM_GHS(arma::mat X, arma::mat S, arma::mat theta, arma::mat sigma, arma::mat Lambda_sq, double epsilon, bool verbose, int maxitr, bool savepath, bool save_Q, double tau_sq, double machine_eps, bool fix_tau = false, bool stop_underflow=false, bool savepath_tau=false) {

  // get dimensions
  const int M = X.n_cols;
  const int N = X.n_rows;
  
  // For saving variables
  int save_dim;
  if (M < 201 & savepath==true){ // Saving exhausts memory if M>201
    save_dim = maxitr;
  }
  else{
    save_dim = 1;
  }
  arma::cube theta_path(M, M, save_dim);
  arma::vec tau_sq_all(maxitr);
  arma::vec Q_vals(maxitr);
  double Q_val_old= -numeric_limits<double>::max();
  double Q_val_new;
  
  // initialize intermediate values
  int niter,i, count;
  double eps;
  arma::uvec pseq(M);
  for(i = 0; i < M; i++){
    pseq(i) = i;
  }
  
  // Initialise updates
 arma::mat theta_update = theta; // Save previous estimate to assess convergence
 arma::mat E_Nu_mat(M,M);
 arma::mat E_NuInv(M,M);
 arma::mat E_XiInv(M,M); // Used if variables are grouped
 double E_xiInv; // Used if variables are not grouped
 arma::cube theta_sigma_update(M,M,2);
 eps = epsilon + 1;
 niter = 1;
 count = 0;
 List list;
 
  if(savepath){
    theta_path.slice(0) = theta;
  }
  if(savepath_tau){
    tau_sq_all(0) = tau_sq;
  }
    
  while(eps>epsilon & count < maxitr){
      
    // E-step
    E_NuInv = E_Nu(Lambda_sq);
    E_xiInv = E_xi(tau_sq);

    // M-step
    // Update Lambda_sq, tau and theta in the M-step
    Lambda_sq = M_lambda(N, M, theta, E_NuInv, machine_eps, stop_underflow, tau_sq);
    if (fix_tau == false){
        tau_sq = M_tau(M, theta, Lambda_sq, E_xiInv, machine_eps, stop_underflow);
    }
    theta_sigma_update = M_theta(N, M, theta, S, sigma, Lambda_sq, pseq, machine_eps, stop_underflow, tau_sq); // Pass S as dummy argument
    // Get updates
    theta_update = theta_sigma_update.slice(0);
    sigma = theta_sigma_update.slice(1);

    if(savepath){
      theta_path.slice(count+1) = theta_update;
    }
    if(savepath_tau){
      tau_sq_all(count+1)=tau_sq;
    }
    // Evaluate objective function if it is to be saved
    if (save_Q) {
      Q_val_new = Q_val(N, M, theta, S, Lambda_sq, E_NuInv, tau_sq, E_xiInv); 
      Q_vals(count) = Q_val_new;
      // Update estimate
      Q_val_old = Q_val_new;
    } 
    // Find update max diff
    eps = max(max(abs(theta_update - theta)));
    
    theta = theta_update;
    count++;
    if(verbose){
      Rcout << "Itr = " << count << " Max diff = " << eps << endl;
    }
  }
  theta = theta_update;
  if(save_Q){
    Q_vals= Q_vals.head(count);
  }
  // Save results
  list["S"] = S;
  list["theta"] = theta;  
  list["sigma"] = sigma;
  list["X"] = X;
  list["Lambda_sq"] = Lambda_sq;
  list["tau_sq"] = tau_sq;
  if(save_Q){
    list["Q_vals"] = Q_vals;  
  }
  if(savepath){
    list["theta_path"] = theta_path;
  }
  if(savepath_tau){
    list["count"] = count;
    list["tau_sq_path"] = tau_sq_all;
  }
  return list;
}