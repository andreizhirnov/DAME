#include "RcppArmadillo.h"
#include <memory>

// [[Rcpp::depends(RcppArmadillo)]] 

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// parent class
class Effects {
public:
  arma::mat beta;
  Effects(arma::mat beta): beta(beta) {}
  virtual ~Effects() {}
  virtual arma::mat ym(arma::mat mu) {
    return mu;
  } 
  virtual arma::mat dydm(arma::mat mu) {
    arma::mat r;
    r.ones(arma::size(mu));
    return r; 
  } 
  virtual arma::mat d2ydm2(arma::mat mu) {
    arma::mat r;
    r.zeros(arma::size(mu));
    return r; 
  }
  arma::mat mx(arma::mat const& mmat) {
    return beta*mmat.t();
  }
// discrete
  arma::mat ddx(arma::mat const& mmat_step, arma::mat const& mmat) {
    return ym(mx(mmat_step)) - ym(mx(mmat));
  }
// derivatives
  arma::mat dydx(arma::mat const& mmat, arma::mat const& d2mdxdb) {
    arma::mat dmdx = beta*d2mdxdb.t(); 
    arma::mat dydm_top = dydm(mx(mmat));
    if (dmdx.n_cols==1) {
      dydm_top.each_col() %= dmdx.col(0);
    } else {
      dydm_top %= dmdx;
    }
    return dydm_top;
  }
  
// discrete: matrix for delta method
  arma::mat ddx_L(const arma::mat& mmat_step, const arma::mat& mmat){
    arma::mat h_step = mmat_step.t(); // t by n
    arma::mat h = mmat.t(); // t by n
    h_step.each_row() %= dydm(mx(mmat_step)).row(0);
    h.each_row() %= dydm(mx(mmat)).row(0);
    return h_step - h;
};
// derivatives: matrix for delta method
  arma::mat dydx_L(arma::mat const& mmat,
                   arma::mat const& d2mdxdb) {
    arma::mat h = mmat.t(); // t by n
    arma::mat mxv = mx(mmat);  // i by n
    arma::rowvec h1 = d2ydm2(mxv).row(0); // i by n
    arma::mat h2 = beta*d2mdxdb.t(); // i by n or i by 1
    if (h2.n_cols==1) {
      h1 *= h2(0);
    } else {
      h1 %= h2.row(0);
    }
    h.each_row() %= h1; 
    arma::mat g = d2mdxdb.t(); // t by n or t by 1
    if (g.n_cols==1) {
      g *= dydm(mxv).row(0);
    } else {
      g.each_row() %= dydm(mxv).row(0);
    }  
    return h + g;
  }
};

// c('logit'1, 'probit'2, 'cauchit'3, 'cloglog'4, 'log'5, 'sqrt'6, "1/mu^2"7, "inverse"8)
// logit
class Effects_1 : public Effects {
public:
  Effects_1(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override {
    arma::mat r=mu;
    r.transform( [](double x ) { return R::plogis(x, 0.0, 1.0, 1, 0); } );
    return r;
  }
  arma::mat dydm(arma::mat mu) override {
    arma::mat r=ym(mu);
    return r % (1.0-r);
  }
  arma::mat d2ydm2(arma::mat mu) override {
    arma::mat r=ym(mu);
    return r % (1.0-r) % (1.0-2.0*r);
  }
};

// probit
class Effects_2 : public Effects {
public:
  Effects_2(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override {
    arma::mat r=mu;
    r.transform( [](double x ) { return R::pnorm(x, 0.0, 1.0, 1, 0); } );
    return r;
  }
  arma::mat dydm(arma::mat mu) override {
    arma::mat r=mu;
    r.transform( [](double x ) { return R::dnorm(x, 0.0, 1.0, 0); } );
    return r;
  }
  arma::mat d2ydm2(arma::mat mu) override {
    arma::mat r=dydm(mu);
    return r % mu; 
  }
};

// cauchit
class Effects_3 : public Effects {
public:
  Effects_3(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override {
    arma::mat r=mu;
    r.transform( [](double x ) { return R::pcauchy(x, 0.0, 1.0, 1, 0); } );
    return r;
  }
  arma::mat dydm(arma::mat mu) override {
    arma::mat r=mu;
    r.transform( [](double x ) { return R::dcauchy(x, 0.0, 1.0, 0); } );
    return r;
  }
  arma::mat d2ydm2(arma::mat mu) override {
    arma::mat r;
    r.zeros(arma::size(mu));
    return r;
  }
};

// cloglog
class Effects_4 : public Effects {
public:
  Effects_4(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override { 
    return 1.0-arma::exp(-arma::exp(mu));
  }
  arma::mat dydm(arma::mat mu) override { 
    return arma::exp(mu - arma::exp(mu)); 
  }
  arma::mat d2ydm2(arma::mat mu) override {
    return -arma::exp(mu - arma::exp(mu))*(arma::exp(mu)-1.0);
  }
};

// log
class Effects_5 : public Effects {
public:
  Effects_5(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override { 
    return arma::exp(mu);
  }
  arma::mat dydm(arma::mat mu) override { 
    return arma::exp(mu);
  }
  arma::mat d2ydm2(arma::mat mu) override {
    return arma::exp(mu);
  }
};

// sqrt
class Effects_6 : public Effects {
public:
  Effects_6(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override { 
    return arma::pow(mu,2);
  }
  arma::mat dydm(arma::mat mu) override { 
    return 2.0*mu;
  }
  arma::mat d2ydm2(arma::mat mu) override {
    arma::mat r(arma::size(mu), arma::fill::value(2.0)); 
    return r;
  }
};

// "1/mu^2"
class Effects_7 : public Effects {
public:
  Effects_7(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override { 
    return arma::pow(mu,-0.5);
  }
  arma::mat dydm(arma::mat mu) override { 
    return arma::pow(mu, -1.5)/(-2.0);
  }
  arma::mat d2ydm2(arma::mat mu) override {
    return arma::pow(mu, -2.5)*0.75;
  }
};

// inverse
class Effects_8 : public Effects {
public:
  Effects_8(arma::mat b): Effects(b) {}
  arma::mat ym(arma::mat mu) override { 
    return arma::pow(mu, -1);
  }
  arma::mat dydm(arma::mat mu) override { 
    return -arma::pow(mu, -2);
  }
  arma::mat d2ydm2(arma::mat mu) override {
    return 2.0*arma::pow(mu, -3);
  }
};

// set up an object with coefficients
std::unique_ptr<Effects> make_proc_mc(const arma::vec& beta, 
                     const arma::mat& vcov,
                     const unsigned int& link=0,
                     const unsigned int& iter=1000){
  arma::mat coef = beta;
  if (iter>0) {
    coef.insert_cols(coef.n_cols, arma::mvnrnd(beta, vcov, iter));
  }
  std::unique_ptr<Effects> eff;
  switch(link) {
  case 1: eff = make_unique<Effects_1>(coef.t()); break;
  case 2: eff = make_unique<Effects_2>(coef.t()); break;
  case 3: eff = make_unique<Effects_3>(coef.t()); break;
  case 4: eff = make_unique<Effects_4>(coef.t()); break;
  case 5: eff = make_unique<Effects_5>(coef.t()); break;
  case 6: eff = make_unique<Effects_6>(coef.t()); break;
  case 7: eff = make_unique<Effects_7>(coef.t()); break;
  case 8: eff = make_unique<Effects_8>(coef.t()); break;
  default: eff = make_unique<Effects>(coef.t());
  }
  return eff;
}

std::unique_ptr<Effects> make_proc_delta(const arma::vec& beta,  
                                      const unsigned int& link=0){
  arma::mat coef = beta;
  std::unique_ptr<Effects> eff;
  switch(link) {
  case 1: eff = make_unique<Effects_1>(coef.t()); break;
  case 2: eff = make_unique<Effects_2>(coef.t()); break;
  case 3: eff = make_unique<Effects_3>(coef.t()); break;
  case 4: eff = make_unique<Effects_4>(coef.t()); break;
  case 5: eff = make_unique<Effects_5>(coef.t()); break;
  case 6: eff = make_unique<Effects_6>(coef.t()); break;
  case 7: eff = make_unique<Effects_7>(coef.t()); break;
  case 8: eff = make_unique<Effects_8>(coef.t()); break;
  default: eff = make_unique<Effects>(coef.t());
  }
  return eff;
}

// aggregate
arma::mat agg_mc(const arma::mat& batch,
                     const arma::umat& wei_locs,
                     const arma::vec& wei_vals,
                     const Rcpp::NumericVector& probs){

  arma::mat wei(arma::max(wei_locs.row(0)) + 1, arma::max(wei_locs.row(1)) + 1); /* reconstruct the weights matrix */
  wei.elem(arma::sub2ind(arma::size(wei), wei_locs)) = wei_vals;
  arma::mat agg = batch*wei;
  agg.each_row() /= arma::sum(wei, 0); /* divide by total */
// extract the estimate and the quantiles

  arma::rowvec est = agg.row(0);
  arma::rowvec se = arma::stddev(agg.tail_rows(agg.n_rows-1), 0, 0);
  arma::mat qs = arma::quantile(agg.tail_rows(agg.n_rows-1), Rcpp::as<arma::vec>(probs), 0);
  return (arma::join_cols(est, se, qs).t());
}

arma::mat agg_delta(const arma::rowvec& est,
                     const arma::mat& bun, 
                     const arma::mat& vcov, 
                     const arma::umat& wei_locs,
                     const arma::vec& wei_vals,
                     const Rcpp::NumericVector& probs){
  
  arma::mat wei(arma::max(wei_locs.row(0)) + 1, arma::max(wei_locs.row(1)) + 1); /* reconstruct the weights matrix */
  // apply weights  
  wei.elem(arma::sub2ind(arma::size(wei), wei_locs)) = wei_vals;
  arma::rowvec denom = arma::sum(wei, 0);
  arma::rowvec agg_e = est*wei / denom;
  arma::mat agg_b = bun*wei;
  agg_b.each_row() /= denom;
// extract standard errors and sampling quantiles
  arma::vec stdqs = Rcpp::qnorm(probs);  
  arma::vec se = arma::sqrt(arma::diagvec(agg_b.t()*vcov*agg_b));
  arma::mat quantiles = stdqs*se.t();
  quantiles.each_row() += agg_e;
// extract the estimate and the quantiles  
  
  return (arma::join_cols(agg_e, se.t(), quantiles).t()); 
}

// combine the pieces into the calculations of marginal effects

// [[Rcpp::export]]
arma::mat get_ddx_mc(const arma::vec& beta, const arma::mat& vcov,
               const arma::mat& x,
               const arma::mat& xp, 
               const arma::umat& wei_locs,
               const arma::vec& wei_vals,
               const Rcpp::NumericVector& probs,
               const unsigned int& link=0,
               const unsigned int& iter=1000){
  std::unique_ptr<Effects> eff = make_proc_mc(beta, vcov, link, iter);
  arma::mat batch = eff->ddx(xp,x);
  return agg_mc(batch, wei_locs, wei_vals, probs);
}

// [[Rcpp::export]]
arma::mat get_dydx_mc(const arma::vec& beta, const arma::mat& vcov,
                      const arma::mat& x,
                      const arma::mat & d2mdxdb,
                      const arma::umat& wei_locs,
                      const arma::vec& wei_vals,
                      const Rcpp::NumericVector& probs,
                      const unsigned int& link=0,
                      const unsigned int& iter=1000){
  std::unique_ptr<Effects> eff = make_proc_mc(beta, vcov, link, iter); 
  arma::mat batch = eff->dydx(x, d2mdxdb);
  return agg_mc(batch, wei_locs, wei_vals, probs);
}

// [[Rcpp::export]]
arma::mat get_ddx_delta(const arma::vec& beta, const arma::mat& vcov,
                      const arma::mat& x,
                      const arma::mat& xp, 
                      const arma::umat& wei_locs,
                      const arma::vec& wei_vals,
                      const Rcpp::NumericVector& probs,
                      const unsigned int& link=0){
  std::unique_ptr<Effects> eff = make_proc_delta(beta, link);
  arma::rowvec est = eff->ddx(xp,x).row(0);
  arma::mat bun = eff->ddx_L(xp,x);
  return agg_delta(est, bun, vcov, wei_locs, wei_vals, probs);
}

// [[Rcpp::export]]
arma::mat get_dydx_delta(const arma::vec& beta, const arma::mat& vcov,
                         const arma::mat& x,
                         const arma::mat & d2mdxdb,
                         const arma::umat& wei_locs,
                         const arma::vec& wei_vals,
                         const Rcpp::NumericVector& probs,
                         const unsigned int& link=0){
  std::unique_ptr<Effects> eff = make_proc_delta(beta, link);
  arma::rowvec est = eff->dydx(x, d2mdxdb).row(0);
  arma::mat bun = eff->dydx_L(x, d2mdxdb);
  return agg_delta(est, bun, vcov, wei_locs, wei_vals, probs);
}

// [[Rcpp::export]]
arma::vec count_nearest(const arma::mat& vals,
                        const arma::mat& grid,
                        const arma::vec& weights)
{
  const bool weighted = (weights.n_elem == vals.n_rows);
  const double denom = weighted ? arma::mean(weights) : 1.0;
  
  arma::vec wei(grid.n_rows, arma::fill::zeros);
  
  const arma::uword n_grid = grid.n_rows;
  const arma::uword dim = grid.n_cols;
  
  for (arma::uword j = 0; j < vals.n_rows; ++j) {
    const arma::rowvec v = vals.row(j);
    
    arma::uword best = 0;
    double best_dist = std::numeric_limits<double>::infinity();
    for (arma::uword i = 0; i < n_grid; ++i) {
      const arma::rowvec g = grid.row(i);
      double d = arma::vecnorm(g-v);
      if (d < best_dist) {
        best_dist = d;
        best = i;
      }
    } 
    wei(best) += weighted ? weights(j) : 1.0;
  } 
  return wei / denom;
}


