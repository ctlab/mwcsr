#include "include/Parameters.h"

using std::string;
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;
using Rcpp::LogicalVector;
using Rcpp::as;

int getInt(List& params, string name) {
    return as<IntegerVector>(params[name])[0];
}

bool getBool(List& params, string name) {
    return as<LogicalVector>(params[name])[0];
}

Parameters::Parameters(List params) {
    timelimit = getInt(params, "timelimit");
    maxIter = getInt(params, "max_iterations");
    betaIter = getInt(params, "beta_iterations");
    separation = getInt(params, "separation");
    maxAge = getInt(params, "max_age");
    startcons = getBool(params, "start_constraints");
    pegging = getBool(params, "pegging");
    sepIter = getInt(params, "sep_iterations");
    sepIterFreeze = getInt(params, "sep_iter_freeze");
    heurIter = getInt(params, "heur_iterations");
    subgradient = getInt(params, "subgradient");
    beta = as<NumericVector>(params["beta"])[0];
    outputlag = getInt(params, "verbose");
}
