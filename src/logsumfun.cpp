#include <RcppArmadillo.h>
#include "logsumfun.h"

//const double epsilon = std::numeric_limits<double>::min();
const double epsilon = 1e-16;

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
