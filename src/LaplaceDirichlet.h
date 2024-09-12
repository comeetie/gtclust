#ifndef LAPLACEDIRICHLET
#define LAPLACEDIRICHLET



#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;


NumericMatrix outer(const NumericVector & beta){
  int D = beta.length();
  NumericMatrix O(D,D);
  for (int d=0;d<D;d++){
    O(d,_)=beta*beta[d];
  }
  return O;
}


NumericMatrix rownorm(NumericMatrix X,const NumericVector & normfact){
  int D = normfact.length();
  NumericMatrix O(D,D);
  for (int d=0;d<D;d++){
    O(d,_)=X(d,_)/normfact(d);
  }
  return O;
}


NumericMatrix colnorm(NumericMatrix X,const NumericVector & normfact){
  int D = normfact.length();
  NumericMatrix O(D,D);
  for (int d=0;d<D;d++){
    O(_,d)=X(_,d)/normfact(d);
  }
  return O;
}


NumericVector update(const NumericMatrix & H,const NumericVector & g){
  int D = g.length();
  NumericVector v(D);
  for (int d=0;d<D;d++){
    v(d)=sum(H(d,_)*g);
  }
  return v;
}


double lC(const NumericVector & alpha){
  return lgamma(sum(alpha))-sum(lgamma(alpha));
}

//[[Rcpp::export]]
double ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  return sum(exp(beta)*(lpi-lambda))+sum(log(lambda))+n*lC(exp(beta)+1)+sum(beta);
}
//[[Rcpp::export]]
NumericVector grad_ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector alpha = exp(beta);
  NumericVector sa = {sum(alpha+1)};
  double dig_sa =  digamma(sa)[0];
  return alpha*(n*dig_sa-n*Rcpp::digamma(alpha+1)+lpi-lambda)+1;
}

//[[Rcpp::export]]
NumericVector dH_ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector g = grad_ldirpost(beta,n,lambda,lpi);
  return g-n*exp(2*beta)*trigamma(exp(beta)+1)-1;
}

//[[Rcpp::export]]
double z_ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector sa = {sum(exp(beta)+1)};
  double ta =trigamma(sa)[0];
  return n*ta;
}

//[[Rcpp::export]]
double logdetH_ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector dH = dH_ldirpost(beta,n,lambda,lpi);
  double z = z_ldirpost(beta,n,lambda,lpi);
  return log(1+z*sum(exp(2*beta)/dH))+sum(log(abs(dH)));
}

//[[Rcpp::export]]
NumericMatrix iH_ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector dH = dH_ldirpost(beta,n,lambda,lpi);
  double z = z_ldirpost(beta,n,lambda,lpi);
  double nf = 1/z+sum(exp(2*beta)/dH);
  NumericMatrix iH = -1/nf*rownorm(colnorm(outer(exp(beta)),dH),dH);
  int D = beta.length();
  for (int d=0;d<D;d++){
    iH(d,d)=iH(d,d)+1/dH(d);
  }
  return iH;
}

//[[Rcpp::export]]
NumericMatrix H_ldirpost(const NumericVector & beta, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector dH = dH_ldirpost(beta,n,lambda,lpi);
  double z = z_ldirpost(beta,n,lambda,lpi);
  NumericMatrix H = z*outer(exp(beta));
  int D = beta.length();
  for (int d=0;d<D;d++){
    H(d,d)=H(d,d)+dH(d);
  }
  return H;
}

//[[Rcpp::export]]
List ldirpost_norm(const NumericVector beta0, int n,const NumericVector & lambda,const NumericVector & lpi){
  NumericVector beta=clone(beta0);
  int nbmaxit = 50;
  int it = 0;
  int D = beta0.length(); 
  double prec= 1e-6;
  double lp_old;
  NumericVector beta_old;
  double lp=ldirpost(beta,n,lambda,lpi);
  bool conv=FALSE;
  while(!conv){
    lp_old   = lp;
    beta_old = clone(beta);
    NumericMatrix iH=iH_ldirpost(beta,n,lambda,lpi);
    NumericVector g=grad_ldirpost(beta,n,lambda,lpi);
    beta=beta-update(iH,g);
    lp = ldirpost(beta,n,lambda,lpi);
    it = it+1;
    if(Rcpp::traits::is_nan<REALSXP>(lp)){
      Rcout << "Nan value in LL" << std::endl;
      Rcout << beta << std::endl;
      Rcout << lpi << std::endl;
      Rcout << n << std::endl;
      beta=beta_old;
      conv=TRUE;
    }
    if(lp<lp_old && abs(lp-lp_old)>prec){
      beta=beta_old;
      conv=TRUE;
      Rcout << "convergence problem, decrease in L." << std::endl;
      Rcout << abs(lp-lp_old) << std::endl;
      Rcout << beta0 << std::endl;
      Rcout << lpi << std::endl;
      Rcout << n << std::endl;
    }
    if(abs(lp-lp_old) < prec){
      conv=TRUE;
    }
    if(it>nbmaxit){
      conv=TRUE;
      Rcout << "convergence problem, nb maxit reached." << std::endl;
      Rcout << beta0 << std::endl;
      Rcout << lpi << std::endl;
      Rcout << n << std::endl;
    }
  }
  
  double lnf = ldirpost(beta,n,lambda,lpi)+D/2*log(2*M_PI)-0.5*logdetH_ldirpost(beta,n,lambda,lpi);
  List res = List::create(Named("beta",beta),Named("lnf", lnf),Named("nbit",it));
  return res;
}

//[[Rcpp::export]]
double dirichlet_evidence(int n,const NumericVector & lambda,const NumericVector & lpi,const NumericVector & pi){
  int d = pi.length();
  double s_nocstr = n*(d-1)/(-2*sum(pi*((lpi-lambda)-n*log(pi))));
  double s_cstr = 100/min(pi);
  double s = std::max(s_nocstr,s_cstr);
  NumericVector beta0 = log(s*pi-1);
  List res = ldirpost_norm(beta0,n,lambda,lpi);
  return res["lnf"];
}

#endif
