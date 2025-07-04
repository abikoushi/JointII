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

arma::sp_mat rcate(arma::sp_mat & A) {
  std::vector<double> vals;
  std::vector<arma::uword> rows, cols;
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
    vals.push_back(*it);
    rows.push_back(it.row());
    cols.push_back(it.col());
  }
  arma::vec prob(vals);
  arma::uvec row_vec(rows);
  arma::uvec col_vec(cols);
  int K = prob.n_elem;
  int pos;
  arma::vec cumprob = cumsum(prob / sum(prob));
    double u = R::runif(0.0, 1.0);
    for (int k = 0; k < K; ++k) {
      if (u < cumprob(k)) {
        pos = k;
        break;
      }
  }
  arma::sp_mat out(A.n_rows, A.n_cols);
  out(row_vec(pos), col_vec(pos)) = 1.0;
  return out;
}

arma::sp_mat asample(int le_rank, int re_rank, int ls_rank, int rs_rank,
                    const arma::vec& h) {
  int m = h.n_elem;
  arma::sp_mat a(m, m);
  for(int i = le_rank; i < re_rank; ++i){
    int start = std::max(ls_rank, i + 1);
    for(int j = start; j < rs_rank; ++j){
      a(i,j) += h[j-i];
    }
  }
  double den = arma::accu(a);
  if(den > 0.0){
    a /= den;    
  }
  return rcate(a);
}


void Asample(arma::sp_mat & a,
             const arma::uvec& le_rank, const arma::uvec& re_rank,
             const arma::uvec& ls_rank, const arma::uvec& rs_rank,
             const arma::vec& h) {
  int m = h.n_elem;
  arma::sp_mat A(m, m);
  for(int i = 0; i < rs_rank.n_elem; ++i){
    A += asample(le_rank(i), re_rank(i), ls_rank(i), rs_rank(i), h);
  }
  a = A;
}

double par_sample(arma::vec & alpha_h,
                arma::vec & h,
                const arma::sp_mat & A) {
  int m = h.n_rows;
  alpha_h.fill(0.0);
  double lp = 0.0;
  for(int e=0; e<m-1; e++){
    for(int s=(e+1); s<m; s++){
      alpha_h(s-e) += A(e,s);
      if(h(s-e) > 0.0){
        lp += A(e,s)*log(h(s-e));
      }
    }
  }
  h = prob_sample(alpha_h);
  return lp;
}

// [[Rcpp::export]]
List conditional_DIC_gibbs(const arma::uvec & EL,
                        const arma::uvec & ER,
                        const arma::uvec & SL,
                        const arma::uvec & SR,
                        const int & m,
                        const int & iter) {
  arma::vec alpha_h = arma::ones<arma::vec>(m);
  arma::vec h = alpha_h/sum(alpha_h);
  arma::mat hhist(m,iter);
  arma::sp_mat A(m,m);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  lp(0) = -std::numeric_limits<double>::max();
  for(int it = 0; it < iter; it++){
    Asample(A, EL, ER, SL, SR, h);
    lp(it) = par_sample(alpha_h, h, A);
    hhist.col(it) = h;
  }
  return List::create(_["incubation"]=hhist,
                      _["logprob"]=lp);
}

/*
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
*/