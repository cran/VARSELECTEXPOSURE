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
//'
//' Calculates likelihood from observed outcome data and given covariate data/parameters.
//'
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
  double like = dot(Y,eta) - sum(log(1+exp(eta)));

  return(like);
}





//[[Rcpp::export]]
int Sample2(arma::vec groupprob){
  arma::vec cumprob=groupprob;
  int m=0;


  for(m=1;m<groupprob.n_rows;m++){
    cumprob[m]=cumprob[m]+cumprob[m-1];
  }
  //Now we have the vector of cumulative probabilities, let's draw a random unif
  double U=as_scalar(arma::randu(1));

  int Which=0;


  if(U<cumprob[0]){
    Which=0;
  }else{

    for(m=0;m<(groupprob.n_rows-2);m++){
      if( (U>cumprob[m]) && (U<cumprob[m+1]) ){
        Which=m+1;
      }

    }

    if(U>cumprob[groupprob.n_rows-2]){
      Which=groupprob.n_rows -1;
    }


  }


  return(Which);

}


//[[Rcpp::export]]
int Sample1(int G){
  arma::vec groupprob(G);
  int m=0;
  double G1=G;
  groupprob.zeros();
  groupprob = groupprob + 1/G1;


  return(Sample2(groupprob));

}



//[[Rcpp::export]]
double GET_EFFECT2(double B0,      //Intercept estimate
                      double BE,      //BetaE estimate
                      arma::vec BETA, //Beta estimates
                      arma::mat X     //Covariate Matrix

){
  double n = X.n_rows;

  double pEvent = 0;
  double pNoEvent= 0;

  arma::vec pi0store(n);
  arma::vec pi1store(n);

  double ateout = 0;

  double i = 0;
  for(int i = 0; i < n; i++){
    pEvent = B0 + BE + dot(BETA, X.row(i));
    pNoEvent = B0 + dot(BETA, X.row(i));
    if(pNoEvent > 50){
      pi0store[i] = 1;
    }else{
      pi0store[i] = exp(pNoEvent)/(1 + exp(pNoEvent));
    }

    if(pEvent > 50){
      pi1store[i] = 1;
    }else{
      pi1store[i] = exp(pEvent)/(1 + exp(pEvent));
    }



  }


  ateout = sum(pi1store - pi0store)/n;
  return(ateout);
}











//' Obtains posterior samples from an MCMC algorithm to perform variable selection.
//'
//' Performs posterior sampling from an MCMC algorithm to estimate average treatment effect and posterior probability of
//' inclusion of candidate variables.
//'
//' @param Y Binary outcome vector.
//' @param Z Matrix of covariates including binary exposure variable.
//' @param PIN Prior probability of inclusion of candidate variables.
//' @param MAX_COV Maximum number of covariates in desired model.
//' @param SdBeta Prior standard deviation for generating distrubtion of proposal coefficients.
//' @param NUM_REPS Number of MCMC iterations to perform.
//' @return Posterior distributions of the estimated average treatment effect, the estimates of the coefficients, and
//'         their posterior probability of inclusion.
//' @useDynLib VARSELECTEXPOSURE
//' @export
//[[Rcpp::export]]
List MCMC_LOGIT_KEEP(arma::vec Y,
                     arma::mat Z,
                     double PIN,
                     double MAX_COV,
                     double SdBeta,
                     double NUM_REPS){

  arma::mat X = Z;
  double B = NUM_REPS;


  arma::vec BETA(X.n_cols);
  BETA.zeros();
  BETA[0] = 0.1;
  BETA[1] = 0.1;
  arma::vec cbeta(X.n_cols);
  cbeta.zeros();
  cbeta = cbeta+1;
  arma::vec abeta = cbeta;
  arma::vec nbeta = cbeta;



  arma::vec ETA(X.n_cols);
  ETA.zeros();
  ETA[0] = 1;
  ETA[1] = 1;
  double beta0 = 0;
  double cbeta0 = 1;
  double abeta0 = 1;
  double nbeta0 = 1;
  double sig = 1;
  double beta0new = 0;
  arma::vec BETA0STORE(B);
  BETA0STORE.zeros();
  arma::mat BETASTORE(B, BETA.n_elem);
  arma::mat ETASTORE = BETASTORE;
  arma::vec BETANEW = BETA;
  double asigma = 1;
  double nsigma = 1;
  double csigma = 1;

  //add step
  arma::vec MU(X.n_rows);
  MU.zeros();
  double FIRST_DERIV = 0;
  double SECOND_DERIV = 0;
  double MEAN = 0;
  arma::vec sumvec(X.n_rows);


  //ATE/RTE storage
  arma::vec ATESTORE(B);
  arma::vec RTESTORE(B);

  //decision values
  double alpha = 0;
  double U = 0;

  //loop integers
  int b = 0;
  int j = 0;
  int t = 0;


  for(int b = 0; b < B; b++){

    if(b < (B/2)){

      if(b%100 == 0){
        if((abeta0/nbeta0)>.6){
          cbeta0 = cbeta0*2;
        }

        //If acceptance is too small, we need to decrease the variance
        if((abeta0/nbeta0)<.2){
          cbeta0 = cbeta0/2;
        }

        for(int j = 0; j < BETA.n_elem; j++){
          if((abeta[j]/nbeta[j])>.6){
            cbeta[j] = cbeta[j]*2;
          }

          if((abeta[j]/nbeta[j])<.2){
            cbeta[j] = cbeta[j]/2;
          }
        }


        if((asigma/nsigma)>.6){
          csigma = csigma*2;
        }

        if((asigma/nsigma)<.2){
          csigma = csigma/2;
        }

      }



      //Reset counters
      abeta0 = 1;
      nbeta0 = 1;
      asigma = 1;
      nsigma = 1;
      abeta.zeros();
      abeta = abeta+1;
      nbeta.zeros();
      nbeta = nbeta+1;
    }

    beta0new = as_scalar(arma::randn(1))*cbeta0 + beta0;

    U = log(as_scalar(arma::randu(1)));
    alpha = LIKE(Y,X,beta0new,BETA) - LIKE(Y,X,beta0,BETA);
    if(U < alpha){
      beta0 = beta0new;
      abeta0 = abeta0+1;
      BETA0STORE[b] = beta0new;
    }
    nbeta0 = nbeta0+1;

    for(int j = 0; j < BETA.n_elem; j++){
      if(ETA[j] == 1){
        BETANEW = BETA;
        BETANEW[j] = as_scalar(arma::randn(1))*cbeta[j] + BETA[j];
        alpha = LIKE(Y,X,beta0,BETANEW)-LIKE(Y,X,beta0,BETA)-(pow(BETANEW[j],2))/(2*SdBeta)+(pow(BETA[j],2))/(2*SdBeta);

        U = log(as_scalar(arma::randu(1)));

        if(U<alpha){
          abeta[j] = abeta[j]+1;
          BETA = BETANEW;
        }
        nbeta[j] = nbeta[j]+1;
      }
    }

    //Variable Selection
    double IND = Sample1(BETA.n_elem-1)+1;
    if(ETA[IND] == 1){
      //Propose delete
      BETANEW = BETA;
      BETANEW[IND] = 0;
      alpha = LIKE(Y,X,beta0,BETANEW)-LIKE(Y,X,beta0,BETA) + log(1-PIN)-log(PIN)+pow(BETA[IND],2)/(2*SdBeta);
      U = log(as_scalar(arma::randu(1)));

      if(U<alpha){
        BETA = BETANEW;
        ETA[IND] = 0;
      }

    }else{
      if(sum(ETA)-1 < MAX_COV){
        //Propose add
        BETANEW = BETA;
        MU = exp(X*BETANEW + beta0)/(1 + exp(X*BETANEW + beta0));
        FIRST_DERIV = sum(X.col(IND).t()*(Y-MU));
        for(int t = 0; t < X.n_rows; t++){
          sumvec[t] = X.col(IND)[t]*X.col(IND)[t]*MU[t]*(1-MU)[t];
        }
        SECOND_DERIV = -sum(sumvec);
        MEAN = BETANEW[IND] - FIRST_DERIV/SECOND_DERIV;
        BETANEW[IND] = as_scalar(arma::randn(1))*MEAN + 0.1;
        U = log(as_scalar(arma::randu(1)));


        alpha = LIKE(Y,X,beta0,BETANEW)-LIKE(Y,X,beta0,BETA) + log(PIN)-log(1-PIN)-(pow(BETA[IND],2))/(2*SdBeta);

        if(U<alpha){
          BETA = BETANEW;
          ETA[IND] = 1;
        }
      }
    }



    ATESTORE[b] = GET_EFFECT2(beta0, BETA[0], BETA.subvec(1, BETA.n_elem-1), X.cols(1, X.n_cols-1));
    BETA0STORE[b] = beta0;
    BETASTORE.row(b) = BETA.t();
    ETASTORE.row(b) = ETA.t();
  }



  List out = List::create(ATESTORE.subvec(B/2, B-1), BETA0STORE.subvec(B/2, B-1),
                          BETASTORE.rows(B/2, B-1), ETASTORE.rows(B/2, B-1));
  return(out);
}





