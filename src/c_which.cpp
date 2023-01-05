# include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <numeric>

using namespace Rcpp;
using namespace sugar;
using namespace arma;

arma::vec  c_which(int input_id, arma::vec id,  arma::vec n_obs, NumericVector indata){

  int input_id_R = input_id - 1;
  arma:: vec id_R = id - 1;
  arma:: vec out_index(n_obs[input_id_R]);  out_index.fill(0);
  arma:: vec out_values(n_obs[input_id_R]);  out_values.fill(0);

  arma:: vec index;
  int n_tot = sum(n_obs);

  // first find the index of the order
  int j = 0;
  for (int i = 0; i < n_tot; ++i) {
    if (id[i] == input_id){ out_index[j] = i;
      j = j+1;}}
  // find the indexed values;
  for (int i = 0; i < n_obs[input_id_R]; ++i) {
    out_values[i] = indata[out_index[i]];}
  return out_values;
}
