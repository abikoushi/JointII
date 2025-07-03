// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "logsumfun.h"
using namespace Rcpp;
using namespace arma;

// Utility functions
arma::vec prob_sample(const arma::vec& x) {
  int n = x.n_elem;
  arma::vec out(n);
  for(int i = 0; i < n; ++i){
    out[i] = (x[i] > 0.0) ? R::rgamma(x[i], 1.0) : 0.0;
  }
  double total = arma::sum(out);
  if (total > 0.0)
    out /= total;
  return out;
}

arma::vec freq_sample(const arma::vec& x) {
  int n = x.n_elem;
  arma::vec hs(n);
  for(int i = 0; i < n; ++i){
    hs[i] = (x[i] > 0.0) ? R::rgamma(x[i], 1.0) : 0.0;
  }
  return arma::cumsum(hs);
}


arma::vec prob2ccdf_cpp(const arma::vec& x) {
  return arma::reverse(arma::cumsum(arma::reverse(x)));
}

arma::vec colmarginal_cpp(const arma::mat & x) {
  return arma::cumsum(arma::sum(x, 1));
}

//////////
/// Core sampler logic
//////////

arma::vec rcate(int size, arma::vec prob) {
  int K = prob.n_elem;
  arma::vec out(K, fill::zeros);
  arma::vec cumprob = cumsum(prob / sum(prob));
  for (int i = 0; i < size; ++i) {
    double u = R::runif(0.0, 1.0);
    for (int k = 0; k < K; ++k) {
      if (u < cumprob[k]) {
        out(k) = 1.0;
        break;
      }
    }
  }
  return out;
}

arma::vec asample(int le_rank, int re_rank, int ls_rank, int rs_rank,
                         const arma::vec& h, const arma::vec& lam) {
  int m = lam.n_elem;
  arma::sp_mat a(m, m);
  for(int i = le_rank; i < re_rank; ++i){
    int start = std::max(ls_rank, i + 1);
    for(int j = start; j < rs_rank; ++j){
      a(i,j) += h[j-i] * lam[i];
    }
  }
  double den = arma::accu(a);
  if(den > 0.0){
    a /= den;    
  }
  return ;
}


void Asample(arma::sp_mat& a,
                 const arma::ivec& le_rank, const arma::ivec& re_rank,
                 const arma::ivec& ls_rank, const arma::ivec& rs_rank,
                 const arma::vec& h, const arma::vec& lam) {
  int m = lam.n_elem;
  arma::sp_mat A(m, m);
  for(int i = 0; i < rs_rank.n_elem; ++i){
    A += asample(le_rank[i], re_rank[i], ls_rank[i], rs_rank[i], h, lam);
  }
  a = A;
}


arma::sp_mat bsample(int le_rank, int re_rank, int ls_rank, int rs_rank,
                         const arma::vec& h, const arma::vec& lam) {
  int m = lam.n_elem;
  arma::sp_mat B(m, m);
  arma::vec rho(2, fill::zeros);
  for(int i = le_rank; i <= re_rank; ++i){
    int start = std::max(ls_rank, i + 1);
    for(int j = start; j <= rs_rank; ++j){
      int cutoff = m - (rs_rank - (j - i));
      if (cutoff > 0)
        rho[0] += lam[i] * arma::sum(h.subvec(0, cutoff - 1));
      if (cutoff < m)
        rho[1] += lam[i] * arma::sum(h.subvec(cutoff, m - 1));
    }
  }
  rho /= arma::sum(rho);
  for(int i = le_rank; i <= re_rank; ++i){
    for(int j = rs_rank + 1; j < m; ++j){
      B(i,j) += h[j - i];
    }
  }
  double denb = arma::accu(B);
  if(denb > 0){
    int b = 0;
    if(rho[1] > 0 && rho[1] < 1)
    b = R::rgeom(rho[1]);
    B /= accu(B);
  }
  return B;
}

void Bsample(arma::sp_mat& b,
             const arma::ivec& le_rank, const arma::ivec& re_rank,
             const arma::ivec& ls_rank, const arma::ivec& rs_rank,
             const arma::vec& h, const arma::vec& lam) {
  int m = h.n_elem;
  arma::sp_mat B(m, m);
  for(int i = 0; i < re_rank.n_elem; ++i){
    B += bsample(le_rank[i], re_rank[i], ls_rank[i], rs_rank[i], h, lam);
  }
  b = B;
}

