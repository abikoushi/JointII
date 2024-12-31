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

void Aup(arma::mat & A, const arma::vec & h, const arma::vec & lam,
          const arma::uvec & lerank, const arma::uvec & rerank,
          const arma::uvec & lsrank, const arma::uvec & rsrank) {
  int n = lsrank.n_rows;
  A.fill(0.0);
  for(int i=0; i<n; i++){
    A += a_upsub(h, lam, lerank(i), rerank(i), lsrank(i), rsrank(i));
  }
}

/*
void Bup(){
  
}
 
 function bup(le_rank,re_rank,ls_rank,rs_rank,h,lam)
 m = size(lam,1)
 B = UpperTriangular(zeros(m,m))
 rho = zeros(2)
#convolution
 for i in le_rank:re_rank
 start = max(ls_rank,i+1)
 for j in start:rs_rank
 rho[1] += lam[i]*sum(h[1:(m-(rs_rank-(j-i)))])
 rho[2] += lam[i]*sum(h[(m-(rs_rank-(j-i))+1):end])
 end
 end
 rho ./= sum(rho)
#expectation
 for i in le_rank:re_rank
 for j in (rs_rank+1):m
 B[i,j] += h[j-i]
 end
 end
 sb = sum(B)
 if sb > 0
 B ./= sb
 B .*= (rho[2]/rho[1])
 end
 return B
 end
 
 function Bup!(b,le_rank,re_rank,ls_rank,rs_rank,h,lam)
 m = size(h,1)
 B = zero(b)
 for i in eachindex(re_rank)
 B .+= bup(le_rank[i],re_rank[i],ls_rank[i],rs_rank[i],h,lam)
 end
 copy!(b,B)
 end
 
 function summingup!(α, β, A, B, h, lam)
 lp = 0.
 m = size(A,1)
 alpha = zero(α)
 beta = zero(β)
 for i in 1:(m-1)
 for j in (i+1):m
 alpha[j-i] += A[i,j] + B[i,j]
 beta[i] += A[i,j] + B[i,j]
 if h[j-i]>0
 lp += xlogy(A[i,j]+B[i,j],h[j-i])
 end
 if lam[i]>0
 lp += xlogy(A[i,j]+B[i,j],lam[i])
 end
 end
 end
 copy!(α,alpha)
 copy!(β,beta)
 return lp
 end
 
*/

double par_up(arma::vec & alpha_h, arma::vec & beta_l, arma::vec & h, arma::vec & lam, const arma::mat & A) {
  int m = h.n_rows;
  alpha_h.fill(0);
  beta_l.fill(0);
  for(int e=0; e<m-1; e++){
    //Rprintf("e %d\n", e);
    for(int s=(e+1); s<m; s++){
      //Rprintf("s %d\n", s);
      alpha_h(s-e) += A(e,s);
      beta_l(e) += A(e,s);
    }
  }
  h = alpha_h/sum(alpha_h);
  lam = beta_l/sum(beta_l);
  return sumxlogy(alpha_h, h) + sumxlogy(beta_l, lam);
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
  arma::mat A = arma::zeros<arma::mat>(m,m);
  arma::vec lp = arma::zeros<arma::vec>(maxit+1);
  arma::vec h = alpha_h/sum(alpha_h);
  arma::vec lam = beta_l/sum(beta_l);
  for(int it=1; it<maxit+1; it++){
    Aup(A, h, lam, EL, ER, SL, SR);
    lp(it) = par_up(alpha_h, beta_l, h, lam, A);
    if(abs(lp(it)-lp(it-1)) < tol){
      lp = lp.rows(1,it);
      break;
    }
  }
  return List::create(_["alpha"]=alpha_h,
                      _["beta"]=beta_l,
                      _["event"]=A,
                      _["lp"]=lp);
}
