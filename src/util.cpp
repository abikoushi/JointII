#include <RcppArmadillo.h>
#include "util.h"

double sumxlogy(const arma::vec & x, const arma::vec & y){
  double out=0;
  for(int i=0; i<x.n_rows; i++){
    if(arma::as_scalar(x.row(i))>0){
      out += arma::as_scalar(x.row(i)*log(y.row(i)));
    }
  }
  return out;
}

arma::vec expmeanlog(arma::vec x){
  double sumx = sum(x);
  arma::vec out = arma::zeros<arma::vec>(x.n_rows);
  //out.fill(0);
  for(int i=0; i<x.n_rows; i++){
    if(x(i)>0){
      out.row(i) = exp(R::digamma(x(i)) - R::digamma(sumx)); 
    }
  }
  return out;
}