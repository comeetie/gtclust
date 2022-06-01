#ifndef GTMETHOD
#define GTMETHOD

#include <Rcpp.h>
#include "node.h"
#include "dist.h"

using namespace Rcpp;


namespace GTMethod{

  class GTMethod {
  
  public:
    virtual void init(const NumericMatrix& X)  {};
    virtual node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      return cnode;
    };
    virtual double dist(node node_g,node node_h) {
      return dist_euclidean_squared(node_g.x,node_h.x);
    };
    virtual node merge(int new_id,node node_g,node node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = (node_g.x + node_h.x)/2;
      new_node.size= node_h.size+node_g.size;
      new_node.height = height;
      return new_node;
    };
    virtual ~GTMethod() {};
  };
  
  // WARD
  class ward : public GTMethod {
  public:
    void init(const NumericMatrix& X) {};
    double dist(node node_g,node node_h) {
      double w = static_cast< double >(node_g.size*node_h.size)/static_cast< double >(node_g.size+node_h.size);
      return w*dist_euclidean_squared(node_g.x,node_h.x);
    };
    node merge(int new_id,node node_g,node node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = (node_g.size*node_g.x + node_h.size*node_h.x)/(node_g.size+node_h.size);
      new_node.size= node_h.size+node_g.size;
      new_node.height = height;
      return new_node;
    }
    ward() : GTMethod() {};
  };
  
  // CENTROID
  class centroid : public GTMethod {
  public:
    void init(const NumericMatrix& X) {};
    node merge(int new_id,node node_g,node node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = (node_g.size*node_g.x + node_h.size*node_h.x)/(node_g.size+node_h.size);
      new_node.size= node_h.size+node_g.size;
      new_node.height = height;
      return new_node;
    }
    centroid() : GTMethod() {};
  };
  
  class median : public GTMethod {
  public:
    void init(const NumericMatrix& X) {};
    median() : GTMethod() {};
  };
  
  class chisq : public GTMethod {
  public:
    void init(const NumericMatrix& X) {
      int D = X.ncol();
      int T = sum(X);
      NumericVector wt(D);
      for(int d=0; d<D; ++d){
        wt(d)=T/sum(X(_,d));
      }
      w=wt;
    };
    double dist(node node_g,node node_h) {
      return dist_chisq(node_g.x,node_h.x,w);
    };
    node merge(int new_id,node node_g,node node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = node_g.x + node_h.x;
      new_node.size= node_h.size+node_g.size;
      new_node.height = height;
      return new_node;
    }
    chisq() : GTMethod() {};
  private:
    NumericVector w;
  };
  
  
  
  // BAYES Mixture of Multinomials
  class bayes_mom : public GTMethod {
  public:
    void init(const NumericMatrix& X) {};
    node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      double ldm = log_dirichlet_multinom(x,beta);
      List cstats = List::create(Named("d",0),
                                 Named("pi",0),
                                 Named("L",ldm),
                                 Named("Lt",ldm),
                                 Named("r",0));
      cnode.stats = cstats;
      return cnode;
    };
    double dist(node node_g,node node_h) {
      // avoid symmetry problems due to numerical problems
      if(node_g.id>node_h.id){
        node nt;
        nt = node_g;
        node_g=node_h;
        node_h=nt;
      }
      double dg = node_g.stats["d"];
      double dh = node_h.stats["d"];
      double Ltg = node_g.stats["Lt"];
      double Lth = node_h.stats["Lt"];
      double Lg = node_g.stats["L"];
      double Lh = node_h.stats["L"];
      double dn   = ladd(lgamma(node_g.size+node_h.size),dg+dh);
      double pi   = lgamma(node_g.size+node_h.size)-dn;
      NumericVector x   = node_g.x + node_h.x;
      double L    = log_dirichlet_multinom(x,beta);
      double Lt   = ladd(pi+L,-pi+Ltg+Lth);
      // reverse for min heap
      //return Lt-pi-L;
      return Lg+Lh-L;
    };
    node merge(int new_id,node node_g,node node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = node_g.x + node_h.x;
      new_node.size= node_h.size+node_g.size;
      new_node.height = height;
      
      double dg = node_g.stats["d"];
      double dh = node_h.stats["d"];
      double Ltg = node_g.stats["Lt"];
      double Lth = node_h.stats["Lt"];
      double Lg = node_g.stats["L"];
      double Lh = node_h.stats["L"];
      double d   = ladd(lgamma(node_g.size+node_h.size),dg+dh);
      double pi  = lgamma(node_g.size+node_h.size)-d;
      double L   = log_dirichlet_multinom(new_node.x,beta);
      double Lt  = ladd(pi+L,-pi+Ltg+Lth);
      //double r   = pi+L-Lt;
      double r = L-Lg-Lh;
      List cstats = List::create(Named("d",d),
                                 Named("pi",pi),
                                 Named("L",L),
                                 Named("Lt",Lt),
                                 Named("r",r));
      new_node.stats = cstats;
      return new_node;
    }
    bayes_mom(double beta_val) {
      beta = beta_val;
    };
  private:
    double beta;
  };
  
  
  
  // BAYES Diagonal Mixture Models
  class bayes_dgmm : public GTMethod {
  public:
    void init(const NumericMatrix& X) {

      if(Rcpp::traits::is_nan<REALSXP>(beta)){
        beta = 0.1*sum(colvar(X)); 
      }

      if(Rcpp::traits::is_nan<REALSXP>(mu[0])){
        mu = colmeans(X);
      }

    };
    node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      NumericVector S(x.length());
      double ldm = log_gauss_evidence(x,S,1,kappa,tau,beta,mu);
      List cstats = List::create(Named("d",0),
                                 Named("pi",0),
                                 Named("L",ldm),
                                 Named("Lt",ldm),
                                 Named("r",0),
                                 Named("S",S));
      cnode.stats = cstats;
      return cnode;
    };
    double dist(node node_g,node node_h) {
      // avoid symmetry problems due to numerical problems
      if(node_g.id>node_h.id){
        node nt;
        nt = node_g;
        node_g=node_h;
        node_h=nt;
      }
      double dg = node_g.stats["d"];
      double dh = node_h.stats["d"];
      double Ltg = node_g.stats["Lt"];
      double Lth = node_h.stats["Lt"];
      
      double Lg = node_g.stats["L"];
      double Lh = node_h.stats["L"];
      double dn   = ladd(lgamma(node_g.size+node_h.size),dg+dh);
      double pi   = lgamma(node_g.size+node_h.size)-dn;
      int size = node_g.size+node_h.size;
      NumericVector x   = node_g.x*node_g.size/size + node_h.x*node_h.size/size;
      NumericVector Sg = node_g.stats["S"];
      NumericVector Sh = node_h.stats["S"];
      NumericVector S  = Sg+node_g.size*(node_g.x-x)*(node_g.x-x)+Sh+node_h.size*(node_h.x-x)*(node_h.x-x);
      double L    = log_gauss_evidence(x,S,size,kappa,tau,beta,mu);
      double Lt   = ladd(pi+L,-pi+Ltg+Lth);
      // reverse for min heap
      //return Lt-pi-L;
      return Lg+Lh-L;
    };
    node merge(int new_id,node node_g,node node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.size= node_h.size+node_g.size;
      new_node.x  = node_g.x*node_g.size/new_node.size + node_h.x*node_h.size/new_node.size;
      new_node.height = height;
      
      double dg = node_g.stats["d"];
      double dh = node_h.stats["d"];
      double Ltg = node_g.stats["Lt"];
      double Lth = node_h.stats["Lt"];
      double Lg = node_g.stats["L"];
      double Lh = node_h.stats["L"];
      NumericVector Sg = node_g.stats["S"];
      NumericVector Sh = node_h.stats["S"];
      NumericVector S  = Sg+node_g.size*(node_g.x-new_node.x)*(node_g.x-new_node.x)+
        Sh+node_h.size*(node_h.x-new_node.x)*(node_h.x-new_node.x);
      
      double d   = ladd(lgamma(node_g.size+node_h.size),dg+dh);
      double pi  = lgamma(node_g.size+node_h.size)-d;
      double L   = log_gauss_evidence(new_node.x,S,new_node.size,kappa,tau,beta,mu);
      double Lt  = ladd(pi+L,-pi+Ltg+Lth);
      double r   = L-Lg-Lh;
      //double r   = pi+L-Lt;
      List cstats = List::create(Named("d",d),
                                 Named("pi",pi),
                                 Named("L",L),
                                 Named("Lt",Lt),
                                 Named("r",r),
                                 Named("S",S));
      new_node.stats = cstats;
      return new_node;
    }
    bayes_dgmm(double kappa_val,double tau_val,double beta_val,NumericVector mu_val) {
      kappa = kappa_val;
      tau   = tau_val;
      beta  = beta_val;
      mu    = mu_val; 
    };
  private:
    double kappa;
    double tau;
    double beta;
    NumericVector mu;
  };

}


#endif
