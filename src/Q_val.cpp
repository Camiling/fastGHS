// [[Rcpp::depends(RcppArmadillo)]]
#include "Q_val.h"
#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

double Q_val(int N, int M, mat &theta, mat&S, mat &Lambda_sq, mat &E_NuInv, double tau_sq, double E_xiInv) {
  double res;
  mat mat_temp;
  double sum_temp;
  mat_temp = -2*log(Lambda_sq) - pow(theta,2)/Lambda_sq/(2*tau_sq) - E_NuInv/Lambda_sq;
  sum_temp = (sum(sum(mat_temp)) - sum(mat_temp.diag()))/2;
  res = M/2*log(det(theta)) - trace(S*theta)/2 + sum_temp - (M*(M-1)/2+3)*log(tau_sq)/2 - E_xiInv/tau_sq;
  return res;
}
