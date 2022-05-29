#ifndef _fastGHS_ECM_GHS_H
#define _fastGHS_ECM_GHS_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


List ECM_GHS(mat X, mat S, mat theta, mat sigma, mat Lambda_sq, double epsilon, bool verbose, int maxitr, bool savepath, bool save_Q, double tau_sq, double machine_eps, bool fix_tau = false, bool stop_underflow=true, bool savepath_tau=false);
#endif