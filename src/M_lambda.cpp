// [[Rcpp::depends(RcppArmadillo)]]
#include "M_lambda.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;


mat M_lambda(int N, int M, mat &theta,mat E_Nu, double machine_eps,bool stop_underflow, double tau_sq) {
  mat Lambda_sq_new;
  Lambda_sq_new = (E_Nu + exp(2*log(abs(theta))- log(2*tau_sq)))/2;
  int i; 
  int j; 
  if(stop_underflow){
    if(min(min(Lambda_sq_new)) < machine_eps){
      for(i=0; i<M; i++){
        for(j=0; j<M; j++){
          if(Lambda_sq_new(i,j) < machine_eps){
            Lambda_sq_new(i,j) = machine_eps;
          }
        }
      }
    } 
  }
  
  return Lambda_sq_new; 
}

