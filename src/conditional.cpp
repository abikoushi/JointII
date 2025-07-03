#include <RcppArmadillo.h>
#include "logsumfun.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::sp_mat a_up0sub(const arma::vec & h,
                     const int & lerank, const int & rerank,
                     const int & lsrank, const int & rsrank) {
  arma::sp_mat a(h.n_rows, h.n_rows);
  for(int e = lerank; e < rerank; e++){
    int sta = std::max(lsrank, e+1);
    for(int s = sta; s < rsrank; s++){
      a(e,s) += h(s-e);
    }
  }
  double den = arma::accu(a);
  if(den > 0.0){
    a /= den;    
  }
  return a;
}

void Aup0(arma::sp_mat & A, 
         const arma::vec & h,
         const arma::uvec & lerank, const arma::uvec & rerank,
         const arma::uvec & lsrank, const arma::uvec & rsrank) {
  A = arma::sp_mat(A.n_rows, A.n_cols);
  for(int i=0; i<lsrank.n_rows; i++){
    A += a_up0sub(h,lerank(i), rerank(i), lsrank(i), rsrank(i));
  }
}

double par_up0(arma::vec & alpha_h,
              arma::vec & h,
              const arma::sp_mat & A) {
  int m = h.n_rows;
  alpha_h.fill(0);
  double lp = 0;
  for(int e=0; e < (m-1); e++){
    for(int s= (e+1); s < m; s++){
      alpha_h(s-e) += A(e,s);
    }
  }
  h = alpha_h/sum(alpha_h);
  lp = sumxlogy(alpha_h, h);
  return lp;
}

// [[Rcpp::export]]
List conditional_DIC_em(const arma::uvec & EL,
                  const arma::uvec & ER,
                  const arma::uvec & SL,
                  const arma::uvec & SR,
                  const int & m,
                  const int & maxit,
                  const double & tol) {
  arma::vec alpha_h = arma::ones<arma::vec>(m);
  arma::vec h = alpha_h/sum(alpha_h);
  h = h/sum(h);
  arma::sp_mat A(m,m);
  arma::vec lp = arma::zeros<arma::vec>(maxit+1);
  lp(0) = -std::numeric_limits<double>::max();
  for(int it=1; it<maxit+1; it++){
    Aup0(A, h, EL, ER, SL, SR);
    lp(it) = par_up0(alpha_h, h, A);
    if(abs(lp(it)-lp(it-1)) < tol){
      lp = lp.rows(0,it);
      break;
    }
  }
  return List::create(_["incubation"]=h,
                      _["alpha"] = alpha_h,
                      _["event"]=A,
                      _["logprob"]=lp);
}

double vpar_up0(arma::vec & alpha_h,
               arma::vec & h,
               const arma::sp_mat & A) {
  int m = h.n_rows;
  alpha_h.fill(0.0);
  double lp = 0.0;
  for(int e=0; e<m-1; e++){
    //Rprintf("e %d\n", e);
    for(int s=(e+1); s<m; s++){
      //Rprintf("s %d\n", s);
      alpha_h(s-e) += A(e,s);
      if(h(s-e) > 0.0){
        lp += A(e,s)*log(h(s-e));
      }
    }
  }
  h = expmeanlog(alpha_h);
  //logh = meanlog(alpha_h);
  //loglam = meanlog(alpha_l);
  //lp = sumxy(alpha_h, logh) + sumxy(alpha_l, loglam);
  //double lp = sumxlogy(alpha_h, h) + sumxlogy(alpha_l, lam);
  return lp;
}

// [[Rcpp::export]]
List conditional_DIC_vb(const arma::uvec & EL,
                  const arma::uvec & ER,
                  const arma::uvec & SL,
                  const arma::uvec & SR,
                  const int & m,
                  const int & maxit,
                  const double & tol) {
  arma::vec alpha_h = arma::ones<arma::vec>(m);
  arma::vec h = alpha_h/sum(alpha_h);
  arma::sp_mat A(m,m);
  arma::vec lp = arma::zeros<arma::vec>(maxit+1);
  lp(0) = -std::numeric_limits<double>::max();
  for(int it=1; it<maxit+1; it++){
    Aup0(A, h, EL, ER, SL, SR);
    lp(it) = vpar_up0(alpha_h, h, A);
    if(abs(lp(it)-lp(it-1)) < tol){
      lp = lp.rows(0,it);
      break;
    }
  }
  return List::create(_["incubation"]=alpha_h/sum(alpha_h),
                      _["alpha"]=alpha_h,
                      _["event"]=A,
                      _["logprob"]=lp);
}
