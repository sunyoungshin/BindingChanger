#ifndef __MOTIF_SCORE__
#define __MOTIF_SCORE__

#include "helper.h"

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
using namespace Rcpp;

RcppExport SEXP motif_score(SEXP, SEXP);
int get_min_start_pos(int, int);
int get_max_start_pos(int, int);
IntegerVector revert_sequence(IntegerVector);
NumericVector comp_subseq_scores(NumericMatrix, IntegerVector);
SequenceScores comp_seq_scores(NumericMatrix, IntegerVector, bool);
double pwm_log_prob(NumericMatrix, IntegerVector, int);
double bidir_pwm_log_prob(NumericMatrix, IntegerVector, int);
RcppExport SEXP transition_matrix(SEXP);
#endif
