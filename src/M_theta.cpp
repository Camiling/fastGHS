#include "M_theta.h"
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h> 

using namespace std;
using namespace arma;
using namespace Rcpp;

cube M_theta(int N, int M, mat theta, mat &S, mat sigma, mat &Lambda_sq, uvec pseq, double machine_eps,bool stop_underflow, double tau_sq) {
  
  // Return a MxMx2 cube with theta and sigma
  cube res(M,M,2);
  mat theta_new = theta;
  mat sigma_new = sigma;
  int i; int j;
  uvec remove_i(M-1); // All indices except i
  uvec left_i(1); // The removed index
  mat theta_mi_mi(M-1, M-1); // Partioned precision matrix
  mat theta_i_mi(M-1, 1); // Vector from partioned precision matrix
  mat theta_i_i(1, 1); // Vector from partioned precision matrix
  mat sigma_mi_mi(M-1, M-1); // Partioned covariance matrix
  mat sigma_i_mi(M-1, 1); // Vector from partioned covariance matrix
  mat sigma_i_i(1,1); // Partioned covariance matrix
  mat S_i_mi(M-1, 1); // Vector from partioned sample matrix
  mat S_i_i(1,1); // element from partioned sample matrix
  mat Lambda_sq_i_mi(M-1, 1); // Vector from partioned Lambda_sq matrix
  mat Tau_sq_i_mi(M-1,1); // Vector of tau_sq values of partition. Not used if varibles are not grouped. 
  
  // Diagonal matrices
  mat Lambda_diag  = zeros<mat>(M-1,M-1);
  mat Tau_diag = zeros<mat>(M-1,M-1);
  
  // Some additional quantities
  mat theta_prod_val(1,1); // Double check dim
  mat theta_prod_vec(M-1,1);
  
  mat theta_mi_mi_Inv(M-1,M-1);
  for(i = 0; i < M; i++){
    remove_i = find(pseq != i);
    left_i(0) = i;
    theta_mi_mi = theta.submat(remove_i, remove_i);
    theta_i_mi = theta.submat(remove_i, left_i);
    sigma_mi_mi = sigma.submat(remove_i, remove_i);
    sigma_i_mi = sigma.submat(remove_i, left_i);
    sigma_i_i = sigma.submat(left_i, left_i);
    S_i_mi = S.submat(remove_i, left_i); 
    S_i_i = S.submat(left_i, left_i); 
    Lambda_sq_i_mi = Lambda_sq.submat(remove_i,left_i);
    // Find inverse of theta submatrix
    theta_mi_mi_Inv = sigma_mi_mi - sigma_i_mi*sigma_i_mi.t()/sigma_i_i(0,0);
    
    Lambda_diag.diag() = 1/Lambda_sq_i_mi;
    theta_i_mi = -inv(S_i_i(0,0)*theta_mi_mi_Inv + Lambda_diag/tau_sq)*S_i_mi;
    theta_i_i = theta_i_mi.t()*theta_mi_mi_Inv*theta_i_mi + M/S_i_i(0,0);
     
    // Avoid computing these quantities multiple times
    theta_prod_vec = theta_mi_mi_Inv*theta_i_mi;
    theta_prod_val = theta_i_i - theta_i_mi.t()*theta_prod_vec;

    if(stop_underflow){
      // Check that no elements are smaller than the machine precision
      if(abs(theta_prod_val(0,0))<machine_eps){
        theta_prod_val(0,0) = sign(theta_prod_val(0,0))*machine_eps;
      }
      for(j=0; j<(M-1); j++){
        if(abs(theta_prod_vec(j, 0))<machine_eps){
          theta_prod_vec(j, 0) = sign(theta_i_mi(j, 0))*machine_eps;
        }
      }
    }
    // Save new sigma  
    sigma.submat(remove_i, remove_i) = theta_mi_mi_Inv + theta_prod_vec*theta_prod_vec.t()/theta_prod_val(0,0);
    sigma.submat(remove_i, left_i) = - theta_prod_vec/theta_prod_val(0,0);
    sigma.submat(left_i,remove_i) = sigma.submat(remove_i, left_i).t();
    sigma.submat(left_i, left_i) = 1/theta_prod_val;
    
    if(stop_underflow){
      // Check that no elements are smaller than the machine precision
      if(abs(theta_i_i(0,0)) < machine_eps){
        theta_i_i(0,0) = sign(theta_i_i(0,0))*machine_eps;
      }
      for(j=0; j<(M-1); j++){
        if(abs(theta_i_mi(j, 0))<machine_eps){
          theta_i_mi(j, 0) = sign(theta_i_mi(j, 0))*machine_eps;
        }
      }
    }
    
    // Save new theta
    theta.submat(remove_i, left_i) = theta_i_mi;
    theta.submat(left_i,remove_i) = theta_i_mi.t();
    theta.submat(left_i,left_i) = theta_i_i;
    
  }
  res.slice(0) = theta;
  res.slice(1) = sigma;
  return(res);
}