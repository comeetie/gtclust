#ifndef GTMETHOD
#define GTMETHOD

#include <Rcpp.h>
#include "node.h"
#include "dist.h"
#include "LaplaceDirichlet.h"

using namespace Rcpp;


namespace GTMethod{

  class GTMethod {
  
  public:
    virtual void init(const NumericMatrix& X)  {};
    virtual abstract_node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      cnode.height = 0;
      return cnode;
    };
    virtual double dist(abstract_node * node_g,abstract_node * node_h) {
      return dist_euclidean_squared(node_g->x,node_h->x);
    };
    virtual node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = (node_g->x + node_h->x)/2;
      new_node.size= node_h->size+node_g->size;
      new_node.height = height;
      return new_node;
    };
    virtual ~GTMethod() {};
  };
  
  // WARD
  class ward : public GTMethod {
  public:
    void init(const NumericMatrix& X) {};
    double dist(abstract_node * node_g,abstract_node * node_h) {
      double w = static_cast< double >(node_g->size*node_h->size)/static_cast< double >(node_g->size+node_h->size);
      return w*dist_euclidean_squared(node_g->x,node_h->x);
    };
    node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = (node_g->size*node_g->x + node_h->size*node_h->x)/(node_g->size+node_h->size);
      new_node.size= node_h->size+node_g->size;
      new_node.height = height;
      return new_node;
    }
    ward() : GTMethod() {};
  };
  
  // CENTROID
  class centroid : public GTMethod {
  public:
    void init(const NumericMatrix& X) {};
    node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = (node_g->size*node_g->x + node_h->size*node_h->x)/(node_g->size+node_h->size);
      new_node.size= node_h->size+node_g->size;
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
    double dist(abstract_node * node_g,abstract_node * node_h) {
      return dist_chisq(node_g->x,node_h->x,w);
    };
    node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = node_g->x + node_h->x;
      new_node.size= node_h->size+node_g->size;
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
    abstract_node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      double ldm = log_dirichlet_multinom(x,beta);
      cnode.height = ldm;
      List cstats = List::create(Named("Lp",ldm));
      cnode.stats = cstats;
      return cnode;
    };
    double dist(abstract_node * node_g,abstract_node * node_h) {
      // just to be sure to avoid symmetry problems due to numerical problems
      if(node_g->id>node_h->id){
        abstract_node * nt;
        nt = node_g;
        node_g=node_h;
        node_h=nt;
      }
      double Lg = node_g->stats["Lp"];
      double Lh = node_h->stats["Lp"];
      NumericVector x   = node_g->x + node_h->x;
      double L    = log_dirichlet_multinom(x,beta)+lgamma(node_g->size+node_h->size+1);
      // reverse since priority queue sorted in ascending order
      return Lg+Lh-L;
    };
    node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = node_g->x + node_h->x;
      new_node.size= node_h->size+node_g->size;
      new_node.height = height;
      double Lp   = log_dirichlet_multinom(new_node.x,beta)+lgamma(new_node.size+1);
      List cstats = List::create(Named("Lp",Lp));
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
        //Rcout << "beta prior was fixed at" << beta << std::endl;
      }

      if(Rcpp::traits::is_nan<REALSXP>(mu[0])){
        mu = colmeans(X);
      }
      //Rcout << "mu prior was fixed at" << mu << std::endl;

    };
    abstract_node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      NumericVector S(x.length());
      double ldm = log_gauss_evidence(x,S,1,kappa,tau,beta,mu);
      List cstats = List::create(Named("Lp",ldm),
                                 Named("S",S));
      cnode.height = ldm;
      cnode.stats = cstats;
      return cnode;
    };
    double dist(abstract_node * node_g,abstract_node * node_h) {
      // avoid symmetry problems due to numerical problems
      if(node_g->id>node_h->id){
        abstract_node * nt;
        nt = node_g;
        node_g=node_h;
        node_h=nt;
      }

      
      double Lg = node_g->stats["Lp"];
      double Lh = node_h->stats["Lp"];
      int size = node_g->size+node_h->size;
      NumericVector x   = node_g->x*node_g->size/size + node_h->x*node_h->size/size;
      NumericVector Sg = node_g->stats["S"];
      NumericVector Sh = node_h->stats["S"];
      NumericVector S  = Sg+node_g->size*(node_g->x-x)*(node_g->x-x)+Sh+node_h->size*(node_h->x-x)*(node_h->x-x);
      double Lp    = log_gauss_evidence(x,S,size,kappa,tau,beta,mu);//+lgamma(size+1);
      // reverse for min heap
      return Lg+Lh-Lp;
    };
    node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.size= node_h->size+node_g->size;
      new_node.x  = node_g->x*node_g->size/new_node.size + node_h->x*node_h->size/new_node.size;
      new_node.height = height;
      NumericVector Sg = node_g->stats["S"];
      NumericVector Sh = node_h->stats["S"];
      NumericVector S  = Sg+node_g->size*(node_g->x-new_node.x)*(node_g->x-new_node.x)+
        Sh+node_h->size*(node_h->x-new_node.x)*(node_h->x-new_node.x);
      double L   = log_gauss_evidence(new_node.x,S,new_node.size,kappa,tau,beta,mu);
      double Lp  = L;//+lgamma(new_node.size+1);
      List cstats = List::create(Named("Lp",Lp),
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

  
  
  // BAYES Dirichlet Mixture Models
  class bayes_dirichlet : public GTMethod {
  public:
    void init(const NumericMatrix& X) {
      if(Rcpp::traits::is_nan<REALSXP>(lambda[0])){
        lambda = colmeans(X);
        lambda.fill(0.01);
      }
    };
    abstract_node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      double ldm = dirichlet_evidence(1,lambda,log(x),x);
      List cstats = List::create(Named("Lp",ldm),Named("lpi",log(x)));
      cnode.stats = cstats;
      cnode.height = ldm;
      return cnode;
    };
    double dist(abstract_node * node_g,abstract_node * node_h) {
      // avoid symmetry problems due to numerical problems
      if(node_g->id>node_h->id){
        abstract_node * nt;
        nt = node_g;
        node_g=node_h;
        node_h=nt;
      }
      
      
      double Lg = node_g->stats["Lp"];
      double Lh = node_h->stats["Lp"];
      int size = node_g->size+node_h->size;
      NumericVector pi   = node_g->x*node_g->size/size + node_h->x*node_h->size/size;
      NumericVector lpig = node_g->stats["lpi"];
      NumericVector lpih = node_h->stats["lpi"];
      NumericVector lpi = lpig+lpih; 
      double Lp    = dirichlet_evidence(size,lambda,lpi,pi)+lgamma(size+1);
      // reverse for min heap
      return Lg+Lh-Lp;
    };
    node merge(int new_id,abstract_node * node_g,abstract_node * node_h,double height) {
      node new_node;
      new_node.id = new_id;
      new_node.size= node_h->size+node_g->size;
      new_node.x  = node_g->x*node_g->size/new_node.size + node_h->x*node_h->size/new_node.size;
      new_node.height = height;
      NumericVector lpig = node_g->stats["lpi"];
      NumericVector lpih = node_h->stats["lpi"];
      NumericVector lpi = lpig+lpih; 
      double L   = dirichlet_evidence(new_node.size,lambda,lpi,new_node.x);
      double Lp  = L+lgamma(new_node.size+1);
      List cstats = List::create(Named("Lp",Lp),
                                 Named("lpi",lpi));
      new_node.stats = cstats;
      return new_node;
    }
    bayes_dirichlet(NumericVector lambda_val) {
      lambda= lambda_val;
    };
  private:
    NumericVector lambda;
  };
  
  
}





#endif
