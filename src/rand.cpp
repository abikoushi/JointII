#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat rand_gamma(const arma::vec & alpha, const int & np) {
  arma::sp_mat out(np, alpha.n_rows);
  for(int i = 0; i < alpha.n_rows; i++){
    if(alpha(i)>0.0){
      out.col(i) = arma::randg(np, arma::distr_param(alpha(i), 1.0)); 
    }
  }
  return out;
}
