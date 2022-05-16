#ifndef DIST
#define DIST

#include <Rcpp.h>
using namespace Rcpp;

float dist_euclidean_squared(NumericVector x1,NumericVector x2){
  return sum((x1-x2)*(x1-x2));
}

float dist_chisq(NumericVector x1,NumericVector x2,NumericVector w){
  if(sum(x1)==0||sum(x2)==0){
    return 0;
  }
  NumericVector x1n = x1 / sum(x1);
  NumericVector x2n = x2 / sum(x2);
  return sum(w*(x1n-x2n)*(x1n-x2n));
}

float log_dirichlet_multinom(NumericVector x,float beta){
  int d = x.length();
  float icl_emiss = lgamma(d*beta)+sum(lgamma(x+beta))-d*lgamma(beta)-lgamma(sum(x+beta));
  return icl_emiss;
}

float ladd(float a, float b){
  float res = 0;
  if(a>b){
    res = a+log(1+exp(b-a));
  }else{
    res = b+log(1+exp(a-b));
  }
  return res;
}



#endif
