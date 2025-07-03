#include <RcppArmadillo.h>
#include "logsumfun.h"

const double epsilon = std::numeric_limits<double>::min();

double sumxlogy(const arma::vec & x, const arma::vec & y){
  double out=0.0;
  for(int i=0; i<x.n_rows; i++){
    //if(x(i)>0.0 && y(i)>0.0){
    if(x(i)>0.0){
      out += x(i)*log(y(i));
    }
  }
  return out;
}

double xlogy(const double & x, const double & y){
  double out=0.0;
    if(x > 0.0){
      out += x*log(y);
    }
  return out;
}

double sumxy(const arma::vec & x, const arma::vec & y){
  double out=0.0;
  for(int i=0; i<x.n_rows; i++){
    if(x(i) > epsilon){
      out += x(i)*y(i);
    }
  }
  return out;
}

arma::vec expmeanlog(arma::vec x){
  double sumx = sum(x);
  arma::vec out = arma::zeros<arma::vec>(x.n_rows);
  for(int i=0; i<x.n_rows; i++){
    if(x(i) > epsilon){
      out(i) = exp(R::digamma(x(i)) - R::digamma(sumx)); 
    }
  }
  return out;
}

arma::vec meanlog(arma::vec x){
  double sumx = sum(x);
  arma::vec out = arma::zeros<arma::vec>(x.n_rows);
  for(int i=0; i<x.n_rows; i++){
    if(x(i) > epsilon){
      out(i) = R::digamma(x(i)) - R::digamma(sumx);
    }
  }
  return out;
}

double logsumexp(arma::vec x){
  double mx = max(x);
  return mx + log(sum(exp(x-mx)));
};





// double summingup_cpp(arma::vec& alpha, arma::vec& beta,
//                      const arma::sp_mat& A, const arma::sp_mat& B,
//                      const arma::vec& h, const arma::vec& lam) {
//   int m = A.n_rows;
//   arma::vec α(m, fill::zeros);
//   arma::vec β(m, fill::zeros);
//   double lp = 0.0;
//   for(int i = 0; i < m-1; ++i){
//     for(int j = i+1; j < m; ++j){
//       double s = A(i,j) + B(i,j);
//       α[j-i] += s;
//       β[i] += s;
//       if(h[j-i] > 0) lp += xlogy(s, h[j-i]);
//       if(lam[i] > 0) lp += xlogy(s, lam[i]);
//     }
//   }
//   alpha = α;
//   beta = β;
//   return lp;
// }
