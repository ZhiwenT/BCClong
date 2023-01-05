#ifndef CWHICH
#define CWHICH

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <numeric>

using namespace Rcpp;
using namespace sugar;
using namespace arma;

arma::vec  c_which(int input_id, arma::vec id,  arma::vec n_obs, NumericVector indata);

#endif
