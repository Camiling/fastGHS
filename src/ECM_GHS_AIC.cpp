//#include "ECM_GHS_AIC.h"
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

#include "E_xi.h"
#include "M_theta.h"
#include "M_tau.h"
#include "E_Nu.h"
#include "M_lambda.h"
#include "AIC.h"
#include "ECM_GHS.h"
#include "Q_val.h"

using namespace std;
using namespace Rcpp;
using namespace arma;


//' Perform GHS with ECM, selecting tau_sq by AIC
//' 
//' This function performs expectation-conditional-maximation or iterated conditional mode estimation for the graphical horseshoe
//' 
//' @param X \eqn{n} by \eqn{p} matrix of data
//' @param S \eqn{p} by \eqn{p} scatter matrix of the data
//' @param theta initial value of the \eqn{p} by \eqn{p} precision matrix
//' @param sigma initial value of the \eqn{p} by \eqn{p} covariance matrix
//' @param Lambda_sq initial value of matrix of squared local shrinkage parameters
//' @param AIC_eps if AIC_selection == TRUE, the convergence tolerance for the AIC convergence assessment
//' @param tau_sq_vec values of tau_sq to consider  
//' @param epsilon tolerance for the convergence assessment
//' @param verbose logical indicator of printing information at each iteration
//' @param maxitr maximum number of iterations
//' @param savepath logical indicator of saving the estimator at each iteration in the ECM algorithm. Only available for p<200
//' @param save_Q should the value of the objective function at each step be saved?
//' @param tau_sq initial value of squared global shrinkage parameter. If exist_group==T, a dummy value should be provided
//' @param machine_eps numerical. The machine precision
//' @param stop_underflow Should tricks to avoid underflow be used?
//' @param stop_overflow Should tricks to avoid overflow be used?
//' 
//' @return A List with resulting ECM estimates, and saved path and objective function convergence information if requested
//' 
//' @export
// [[Rcpp::export]]
List ECM_GHS_AIC(arma::mat X, arma::mat S, arma::mat theta, arma::mat sigma, arma::mat Lambda_sq, double AIC_eps, arma::vec tau_sq_vec, double epsilon, bool verbose, int maxitr, bool savepath, bool save_Q, double tau_sq, double machine_eps, bool stop_underflow=false, bool stop_overflow=false){
  
  // Get dimensions
  const int N = X.n_rows;
  
  // Initialise estimates
  double eps = 1000;
  int count = 0;
  int maxitr_AIC = maxitr;
  if(maxitr>tau_sq_vec.n_elem){
    maxitr_AIC = tau_sq_vec.n_elem;
  }
  arma::vec AIC_scores(maxitr);
  double AIC_prev = 1000000;
  double AIC_update;
  double tau_sq_val;
  List res_GHS;
  List list;
 
  while(eps > AIC_eps & count < maxitr_AIC){
    // Tau value to use
    tau_sq_val = tau_sq_vec(count);
    
    // Use ECM_GHS with current tau_sq value. Set verbose=false and fix_tau=true
    res_GHS = ECM_GHS(X, S, theta, sigma, Lambda_sq, epsilon, false, maxitr, savepath, save_Q, tau_sq_val, machine_eps, true, stop_underflow, stop_overflow);  
    
    // Compute AIC score
    AIC_update = AIC(S, res_GHS["theta"], N, stop_overflow);
    AIC_scores(count) = AIC_update;
    
    // Change in AIC score
    eps=abs(AIC_prev-AIC_update);
    AIC_prev = AIC_update;
    count++;
    
    if(verbose){
      Rcout << "Itr = " << count << " AIC diff = " << eps << endl;
    }
  }
  // Save final model
  list = res_GHS;
  list["ind"] = count; // minus one, but plus one bc of conversion to R
  list["AIC_scores"] = AIC_scores;
  return list;
}