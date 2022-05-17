#ifndef DIST
#define DIST

#include <Rcpp.h>
using namespace Rcpp;

double dist_euclidean_squared(const NumericVector& x1,const NumericVector& x2){
  return sum((x1-x2)*(x1-x2));
}

double dist_chisq(const NumericVector& x1,const NumericVector& x2,const NumericVector& w){
  if(sum(x1)==0||sum(x2)==0){
    return 0;
  }
  NumericVector x1n = x1 / sum(x1);
  NumericVector x2n = x2 / sum(x2);
  return sum(w*(x1n-x2n)*(x1n-x2n));
}

double log_dirichlet_multinom(NumericVector& x,double beta){
  int d = x.length();
  double icl_emiss = lgamma(d*beta)+sum(lgamma(x+beta))-d*lgamma(beta)-lgamma(sum(x+beta));
  return icl_emiss;
}


double log_gauss_evidence(const NumericVector& x,const NumericVector& S,int ng,double kappa,double tau,double beta,const NumericVector mu){
  NumericVector betan = beta +0.5*S + (tau*ng)/(2*(tau+ng))*(x-mu)*(x-mu);
  double taun = tau+ng;
  double kappan = kappa+(ng/2);
  double log_evidence = sum(lgamma(kappan)-lgamma(kappa)+kappa*log(beta)-kappan*log(betan)+0.5*log(tau)-0.5*log(taun)-ng/2*log(2*M_PI));
  return log_evidence;
}


double ladd(double a, double b){
  double res = 0;
  if(a>b){
    res = a+log(1+exp(b-a));
  }else{
    res = b+log(1+exp(a-b));
  }
  return res;
}

NumericVector colmeans(const NumericMatrix& X) {
  int nr = X.nrow(), nc = X.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += X(i, j);
    }
    ans[j] = sum/nr;
  }
  return ans;
}

NumericVector colvar(const NumericMatrix& X) {
  int nr = X.nrow(), nc = X.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    double sumsq = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += X(i, j);
      sumsq += X(i,j)*X(i,j);
    }
    ans[j] = sumsq/nr - (sum/nr)*(sum/nr);
  }
  return ans;
}

#endif
