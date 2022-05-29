#ifndef _fastGHS_ECM_GHS_AIC_H
#define _fastGHS_ECM_GHS_AIC_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


List ECM_GHS_AIC(mat X, mat S, mat theta, mat sigma, mat Lambda_sq, double AIC_eps, vec tau_sq_vec, double epsilon, bool verbose, int maxitr, bool savepath, bool save_Q, double tau_sq, double machine_eps, bool stop_underflow=true, bool stop_overflow=false);
#endif