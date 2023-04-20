// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LIKE
double LIKE(arma::vec Y, arma::mat X, double beta0, arma::vec beta);
RcppExport SEXP _VARSELECTEXPOSURE_LIKE(SEXP YSEXP, SEXP XSEXP, SEXP beta0SEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(LIKE(Y, X, beta0, beta));
    return rcpp_result_gen;
END_RCPP
}
// Sample2
int Sample2(arma::vec groupprob);
RcppExport SEXP _VARSELECTEXPOSURE_Sample2(SEXP groupprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type groupprob(groupprobSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample2(groupprob));
    return rcpp_result_gen;
END_RCPP
}
// Sample1
int Sample1(int G);
RcppExport SEXP _VARSELECTEXPOSURE_Sample1(SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample1(G));
    return rcpp_result_gen;
END_RCPP
}
// GET_EFFECT2
double GET_EFFECT2(double B0, double BE, arma::vec BETA, arma::mat X);
RcppExport SEXP _VARSELECTEXPOSURE_GET_EFFECT2(SEXP B0SEXP, SEXP BESEXP, SEXP BETASEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< double >::type BE(BESEXP);
    Rcpp::traits::input_parameter< arma::vec >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(GET_EFFECT2(B0, BE, BETA, X));
    return rcpp_result_gen;
END_RCPP
}
// MCMC_LOGIT_KEEP
List MCMC_LOGIT_KEEP(arma::vec Y, arma::mat Z, double PIN, double MAX_COV, double SdBeta, double NUM_REPS);
RcppExport SEXP _VARSELECTEXPOSURE_MCMC_LOGIT_KEEP(SEXP YSEXP, SEXP ZSEXP, SEXP PINSEXP, SEXP MAX_COVSEXP, SEXP SdBetaSEXP, SEXP NUM_REPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< double >::type PIN(PINSEXP);
    Rcpp::traits::input_parameter< double >::type MAX_COV(MAX_COVSEXP);
    Rcpp::traits::input_parameter< double >::type SdBeta(SdBetaSEXP);
    Rcpp::traits::input_parameter< double >::type NUM_REPS(NUM_REPSSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC_LOGIT_KEEP(Y, Z, PIN, MAX_COV, SdBeta, NUM_REPS));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_VARSELECTEXPOSURE_LIKE", (DL_FUNC) &_VARSELECTEXPOSURE_LIKE, 4},
    {"_VARSELECTEXPOSURE_Sample2", (DL_FUNC) &_VARSELECTEXPOSURE_Sample2, 1},
    {"_VARSELECTEXPOSURE_Sample1", (DL_FUNC) &_VARSELECTEXPOSURE_Sample1, 1},
    {"_VARSELECTEXPOSURE_GET_EFFECT2", (DL_FUNC) &_VARSELECTEXPOSURE_GET_EFFECT2, 4},
    {"_VARSELECTEXPOSURE_MCMC_LOGIT_KEEP", (DL_FUNC) &_VARSELECTEXPOSURE_MCMC_LOGIT_KEEP, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_VARSELECTEXPOSURE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
