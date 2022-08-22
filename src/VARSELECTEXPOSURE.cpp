#include "RcppArmadillo.h"
#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <iostream>     // std::cout
#include <cmath>
#include <cfloat>
#include <stdio.h>
#include <float.h>
#define PI 3.14159265


// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;




//' Obtains likelihood
//' Calculates likelihood from observed outcome data and given covariate data/parameters.
//' @param Y Binary outcome vector.
//' @param X Matrix of covariates.
//' @param beta0 Intercept parameter.
//' @param beta Vector of covariate parameters of length p.
//' @return Likelihood
//' @useDynLib VARSELECTEXPOSURE
//' @export
//[[Rcpp::export]]
double LIKE(arma::vec Y,   //Vector of binary outcomes
            arma::mat X,   //Matrix of covariates
            double beta0,  //Intercept parameter
            arma::vec beta //Vector of covariate parameters
){
  arma::vec eta = X*beta + beta0;
  arma::vec LOGL = Y*eta - log(1+exp(eta));

  double LIK = sum(LOGL);
  return(LIK);
}







