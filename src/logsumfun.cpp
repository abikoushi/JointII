#include <RcppArmadillo.h>
#include "logsumfun.h"

double sumxlogy(const arma::vec & x, const arma::vec & y){
  double out=0.0;
  for(int i=0; i<x.n_rows; i++){
    if(x(i)>0.0 && y(i)>0.0){
      out += x(i)*log(y(i));
    }
  }
  return out;
}

double sumxy(const arma::vec & x, const arma::vec & y){
  double out=0.0;
  for(int i=0; i<x.n_rows; i++){
    if(x(i)>0.0){
      out += x(i)*y(i);
    }
  }
  return out;
}

arma::vec expmeanlog(arma::vec x){
  double sumx = sum(x);
  arma::vec out = arma::zeros<arma::vec>(x.n_rows);
  //out.fill(0);
  for(int i=0; i<x.n_rows; i++){
    if(x(i)>0.0){
      out(i) = exp(R::digamma(x(i)) - R::digamma(sumx)); 
    }
  }
  return out;
}

arma::vec meanlog(arma::vec x){
  double sumx = sum(x);
  arma::vec out = arma::zeros<arma::vec>(x.n_rows);
  //out.fill(0);
  for(int i=0; i<x.n_rows; i++){
    if(x(i)>0.0){
      out(i) = R::digamma(x(i)) - R::digamma(sumx);
    }
  }
  return out;
}

double logsumexp(arma::vec x){
  double mx = max(x);
  return mx + log(sum(exp(x-mx)));
};
