#ifndef __IS__
#define __IS__

#include "MotifScore.h"

NumericMatrix p_value(
    NumericMatrix,
    NumericVector,
    NumericMatrix,
    NumericVector,
    double,
    int,
    int,
    LoglikType);
double func_delta(NumericMatrix, NumericVector, NumericMatrix, double, int);
double find_theta(NumericMatrix, NumericVector, NumericMatrix, double, int);
IntegerVector importance_sample(NumericMatrix, NumericVector, NumericMatrix, int);
NumericVector compute_sample_score(NumericMatrix, IntegerVector, int, double);
double find_percentile(NumericVector, double);

NumericMatrix gen_utility_matrix(NumericMatrix, NumericMatrix, int, double);
NumericVector gen_prob_start_pos(NumericMatrix, int, NumericVector);

RcppExport SEXP test_find_percentile(SEXP, SEXP);
RcppExport SEXP compute_p_values(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_find_theta(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_func_delta(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_importance_sample(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_compute_sample_score(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_gen_utility_matrix(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_gen_prob_start_pos(SEXP, SEXP, SEXP);

#endif
