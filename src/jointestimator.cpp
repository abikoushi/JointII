#include <RcppArmadillo.h>
#include "logsumfun.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double par_up(arma::vec & alpha_h, arma::vec & beta_l, 
              arma::vec & h, arma::vec & lam, 
              const arma::sp_mat & A, 
              const arma::sp_mat & B) {
  int m = h.n_rows;
  alpha_h.fill(0);
  beta_l.fill(0);
  double lp = 0;
  for(int e=0; e < (m-1); e++){
    for(int s= (e+1); s < m; s++){
      alpha_h(s-e) += A(e,s) + B(e,s);
      beta_l(e) += A(e,s) + B(e,s);
    }
  }
  h = alpha_h/sum(alpha_h);
  lam = beta_l/sum(beta_l);
  lp = sumxlogy(alpha_h, h) + sumxlogy(beta_l, lam);
  return lp;
}

arma::sp_mat a_upsub(const arma::vec & h, const arma::vec & lam,
                     const int & lerank, const int & rerank,
                     const int & lsrank, const int & rsrank) {
  arma::sp_mat a(h.n_rows, h.n_rows);
  for(int e = lerank; e < rerank; e++){
    int sta = std::max(lsrank, e+1);
    for(int s = sta; s < rsrank; s++){
      a(e,s) += h(s-e)*lam(e);
    }
  }
  double den = arma::accu(a);
  if(den > 0.0){
    a /= den;    
  }
  return a;
}

void Aup(arma::sp_mat & A, 
         const arma::vec & h, const arma::vec & lam,
         const arma::uvec & lerank, const arma::uvec & rerank,
         const arma::uvec & lsrank, const arma::uvec & rsrank) {
  A = arma::sp_mat(A.n_rows, A.n_cols);
  for(int i=0; i<lsrank.n_rows; i++){
    A += a_upsub(h, lam, lerank(i), rerank(i), lsrank(i), rsrank(i));
  }
}

arma::sp_mat b_upsub(const arma::vec & h, const arma::vec & lam,
                     const int & lerank, const int & rerank,
                     const int & lsrank, const int & rsrank){
  arma::vec rho = arma::zeros<arma::vec>(2);
  //convolution
  for(int e = lerank; e < rerank; e++){
    int start = std::max(lsrank, e+1);
    for(int s = start; s < rsrank; s++){
      int split = h.n_rows - (rsrank- (s-e));
      rho(0) += lam(e) * sum(h.rows(0, split - 1));
      rho(1) += lam(e) * sum(h.rows(split, h.n_rows - 1));
    }
  }
  //expectation
  arma::sp_mat b = arma::sp_mat(h.n_rows, h.n_rows);
  arma::sp_mat out = arma::sp_mat(h.n_rows, h.n_rows);
  int end = std::min(rsrank, rerank);
  for(int e=lerank; e<end; e++){
    for(int s=(rsrank+1); s<h.n_rows; s++){
      b(e,s) += h(s-e);
    }
  }
  double sb = arma::accu(b);
  double sr = sum(rho);
  if(sb > 0.0){
    b /= sb;
    rho /= sr;
    if(rho(0) > 0.0 && rho(0)< 1.0){
      b *= rho(1)/rho(0);
      out = b;      
    }
  }
  return out;
}

void Bup(arma::sp_mat & B,
         const arma::vec & h, const arma::vec & lam,
         const arma::uvec & lerank, const arma::uvec & rerank,
         const arma::uvec & lsrank, const arma::uvec & rsrank){
  B = arma::sp_mat(B.n_rows, B.n_cols);
  for(int i=0; i<lsrank.n_rows; i++){
    B += b_upsub(h, lam, lerank(i), rerank(i), lsrank(i), rsrank(i));
  }
}

// [[Rcpp::export]]
List joint_DIC_em(const arma::uvec & EL,
                  const arma::uvec & ER,
                  const arma::uvec & SL,
                  const arma::uvec & SR,
                  const int & m,
                  const int & maxit,
                  const double & tol) {
  arma::vec alpha_h = arma::ones<arma::vec>(m);
  arma::vec beta_l = arma::ones<arma::vec>(m);
  //arma::sp_vec alpha_h(m);
  //arma::sp_vec beta_l(m);
  arma::vec h = alpha_h/sum(alpha_h);
  arma::vec lam = beta_l/sum(beta_l);
  h = h/sum(h);
  lam = lam/sum(lam);
  arma::sp_mat A(m,m);
  arma::sp_mat B(m,m); //ghost
  arma::vec lp = arma::zeros<arma::vec>(maxit+1);
  lp(0) = -std::numeric_limits<double>::max();
  for(int it=1; it<maxit+1; it++){
    Aup(A, h, lam, EL, ER, SL, SR);
    //Rprintf("A ");
    Bup(B, h, lam, EL, ER, SL, SR);
    //Rprintf("B ");
    lp(it) = par_up(alpha_h, beta_l, h, lam, A, B);
    if(abs(lp(it)-lp(it-1)) < tol){
      lp = lp.rows(0,it);
      break;
    }
  }
  return List::create(_["incubation"]=h,
                      _["infection"]=lam,
                      _["alpha"] = alpha_h,
                      _["beta"] = beta_l,
                      _["event"]=A,
                      _["truncated_event"]=B,
                      _["logprob"]=lp);
}



double vpar_up(arma::vec & alpha_h, arma::vec & alpha_l, 
               arma::vec & h, arma::vec & lam, 
               const arma::sp_mat & A, const arma::sp_mat & B) {
  int m = h.n_rows;
  alpha_h.fill(0.0);
  alpha_l.fill(0.0);
  double lp = 0.0;
  for(int e=0; e<m-1; e++){
    //Rprintf("e %d\n", e);
    for(int s=(e+1); s<m; s++){
      //Rprintf("s %d\n", s);
      alpha_h(s-e) += A(e,s) + B(e,s);
      alpha_l(e) += A(e,s) + B(e,s);
      if(h(s-e) > 0.0){
        lp += (A(e,s)+B(e,s))*log(h(s-e));
      }
      if(lam(e) > 0.0){
        lp += (A(e,s)+B(e,s))*log(lam(e));
      }
    }
  }
  h = expmeanlog(alpha_h);
  lam = expmeanlog(alpha_l);
  return lp;
}

// [[Rcpp::export]]
List joint_DIC_vb(const arma::uvec & EL,
                  const arma::uvec & ER,
                  const arma::uvec & SL,
                  const arma::uvec & SR,
                  const int & m,
                  const int & maxit,
                  const double & tol) {
  arma::vec alpha_h = arma::ones<arma::vec>(m);
  arma::vec alpha_l = arma::ones<arma::vec>(m);
  arma::vec h = alpha_h/sum(alpha_h);
  arma::vec lam =alpha_l/sum(alpha_l);
  arma::sp_mat A(m,m);
  arma::sp_mat B(m,m); //ghost
  arma::vec lp = arma::zeros<arma::vec>(maxit+1);
  lp(0) = -std::numeric_limits<double>::max();
  for(int it=1; it<maxit+1; it++){
    Aup(A, h, lam, EL, ER, SL, SR);
    //Rprintf("A ");
    Bup(B, h, lam, EL, ER, SL, SR);
    //Rprintf("B ");
    lp(it) = vpar_up(alpha_h, alpha_l, h, lam, A, B);
    if(abs(lp(it)-lp(it-1)) < tol){
      lp = lp.rows(0,it);
      break;
    }
  }
  return List::create(_["incubation"]=alpha_h/sum(alpha_h),
                      _["infection"]=alpha_l/sum(alpha_l),
                      _["alpha"]=alpha_h,
                      _["beta"]=alpha_l,
                      _["event"]=A,
                      _["truncated_event"]=B,
                      _["logprob"]=lp);
}
