double sumxlogy(const arma::vec & x, const arma::vec & y);

double sumxy(const arma::vec & x, const arma::vec & y);

arma::vec expmeanlog(arma::vec x);

arma::vec meanlog(arma::vec x);

double logsumexp(arma::vec x);

double xlogy(const double & x, const double & y);

double par_up(arma::vec & alpha_h, arma::vec & beta_l, 
              arma::vec & h, arma::vec & lam, 
              const arma::sp_mat & A,
              const arma::sp_mat & B);