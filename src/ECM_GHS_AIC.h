#ifndef _fastGHS_ECM_GHS_AIC_H
#define _fastGHS_ECM_GHS_AIC_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


List ECM_GHS_AIC(mat X, mat S, mat theta, mat sigma, mat Lambda_sq, double AIC_eps, vec tau_sq_vec, double epsilon, bool verbose, int maxitr, bool savepath, int exist_group, uvec group, mat N_groups, bool save_Q, double tau_sq, mat Tau_sq, double machine_eps, bool use_ICM = false, bool GHS_like = false, bool stop_underflow=true);
#endif