#ifndef _fastGHS_AIC_H
#define _fastGHS_AIC_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace Rcpp;
using namespace arma;


double AIC(mat S, mat theta, int n);
#endif