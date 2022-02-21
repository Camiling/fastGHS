// [[Rcpp::depends(RcppArmadillo)]]
//#include "AIC.h"
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


//' Compute AIC score of an estimated precision matrix
//' 
//' This function computes the AIC score of an estimated precision matrix
//' 
//' @param S \eqn{p} by \eqn{p} scatter matrix of the data
//' @param theta estimated \eqn{p} by \eqn{p} precision matrix
//' @param n the number of observations
//' @param stop_overflow should overflow be avoided?
//' 
//' @return The AIC score of the estimated precision matrix
//' 
//' @export
// [[Rcpp::export]]
double AIC(arma::mat S, arma::mat theta, int n, bool stop_overflow){
  double AIC_score;
  double d;
  double loglik;
  int p = S.n_cols;
  // Must round the elements of theta
  vec theta_vec = vectorise(theta);
  d = (sum(abs(theta_vec)> 1e-5)-p)/2;
  // The loglikelihood. Must divide scatter matrix by n-1 to get sample covariance matrix.
  if(stop_overflow){
    // Dividing theta by 5 in order to avoid overflow
    loglik = -p*n*log(2*datum::pi)/2 + n*(log(det(theta/5))/2 + log(pow(5,p))/2)+ n*trace(S*(theta/(n-1)))/2;
  }
  else{
    loglik = -p*n*log(2*datum::pi)/2 + n*log(det(theta))/2 + n*trace(S*theta/(n-1))/2;
  }
  AIC_score = -2*loglik + 2*d;
  return AIC_score;
}