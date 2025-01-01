#include <RcppArmadillo.h>
#include "util.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat a_upsub(const arma::vec & h, const arma::vec & lam,
                const int & lerank, const int & rerank,
                const int & lsrank, const int & rsrank) {
  arma::mat a = arma::zeros<arma::mat>(h.n_rows, h.n_rows);
  for(int e=lerank; e<rerank; e++){
    int sta = std::max(lsrank, e+1);
    for(int s=sta; s<rsrank; s++){
      //Rprintf("%d ",s-e);
      a(e,s) += h(s-e)*lam(e);
    }
    //Rprintf("\n");
  }
  a /= arma::accu(a);
  return a;
}

void Aup(arma::sp_mat & A, const arma::vec & h, const arma::vec & lam,
          const arma::uvec & lerank, const arma::uvec & rerank,
          const arma::uvec & lsrank, const arma::uvec & rsrank) {
  int n = lsrank.n_rows;
  //A.fill(0.0);
  A = arma::sp_mat(A.n_rows, A.n_cols);
  for(int i=0; i<n; i++){
    A += a_upsub(h, lam, lerank(i), rerank(i), lsrank(i), rsrank(i));
  }
}


void Bup(arma::sp_mat & B, const arma::vec & h, const arma::vec & lam,
         const arma::uvec & lerank, const arma::uvec & rerank,
         const arma::uvec & lsrank, const arma::uvec & rsrank){
  int N = lsrank.n_rows;
  //B.fill(0.0);
  B = arma::sp_mat(B.n_rows, B.n_cols);
  arma::vec rho = arma::zeros<arma::vec>(2);
  for(int i=0; i<N; i++){
    //Rprintf("%d : ", i);
    //convolution
    for(int e = lerank(i); e<rerank(i); e++){
      int start = lsrank(i);
      start = std::max(start, e+1);
      for(int s=start; s<rsrank(i); s++){
        int split =  h.n_rows-(rsrank(i)-(s-e))-1;
        //Rprintf("%d ", split);
        rho(0) += lam(e)*sum(h.rows(0, split));
        rho(1) += lam(e)*sum(h.rows(split+1, h.n_rows-1));      
      }
      //Rprintf("\n");
  }
  //  Rprintf("Conv ");
  //expectation
  for(int e=lerank(i); e<rerank(i); e++){
  for(int s=(rsrank(i)+1); s<h.n_rows; s++){
  B(e,s) += h(s-e);
  double sb = arma::accu(B);
    if(sb > 0){
      B /= sb;
      B *= rho(1)/rho(0);
    }
    }
  } 
  //Rprintf("Exp ");
  }
}


double par_up(arma::vec & h, arma::vec & lam, const arma::sp_mat & A, const arma::sp_mat & B) {
  int m = h.n_rows;
  h.fill(0);
  lam.fill(0);
  for(int e=0; e<m-1; e++){
    //Rprintf("e %d\n", e);
    for(int s=(e+1); s<m; s++){
      //Rprintf("s %d\n", s);
      h(s-e) += A(e,s) + B(e,s);
      lam(e) += A(e,s) + B(e,s);
    }
  }
  arma::vec h2 = h/sum(h);
  arma::vec lam2 = lam/sum(lam);
  double lp = sumxlogy(h, h2) + sumxlogy(lam, lam2);
  h = h2;
  lam = lam2;
  return lp;
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
  arma::vec h = alpha_h/sum(alpha_h);
  arma::vec lam = beta_l/sum(beta_l);
  arma::sp_mat A(m,m);
  arma::sp_mat B(m,m); //ghost
  arma::vec lp = arma::zeros<arma::vec>(maxit+1);
  lp(0) = -std::numeric_limits<double>::max();
  for(int it=1; it<maxit+1; it++){
    Aup(A, h, lam, EL, ER, SL, SR);
    //Rprintf("A ");
    Bup(B, h, lam, EL, ER, SL, SR);
    //Rprintf("B ");
    lp(it) = par_up(h, lam, A, B);
    if(abs(lp(it)-lp(it-1)) < tol){
      lp = lp.rows(0,it);
      break;
    }
  }
  return List::create(_["incubation"]=h,
                      _["infection"]=lam,
                      _["event"]=A,
                      _["truncated_event"]=B,
                      _["logprob"]=lp);
}
