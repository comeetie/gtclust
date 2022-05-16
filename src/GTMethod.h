#ifndef GTMETHOD
#define GTMETHOD

#include <Rcpp.h>
#include "node.h"
#include "dist.h"

using namespace Rcpp;


namespace GTMethod{

  class GTMethod {
  
  public:
    virtual void init(const NumericMatrix X)  {};
    virtual node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      return cnode;
    };
    virtual float dist(node node_g,node node_h) {
      return dist_euclidean_squared(node_g.x,node_h.x);
    };
    virtual node merge(int new_id,node node_g,node node_h,float height) {
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
    void init(const NumericMatrix X) {};
    float dist(node node_g,node node_h) {
      float w = static_cast< float >(node_g.size*node_h.size)/static_cast< float >(node_g.size+node_h.size);
      return w*dist_euclidean_squared(node_g.x,node_h.x);
    };
    node merge(int new_id,node node_g,node node_h,float height) {
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
    void init(const NumericMatrix X) {};
    node merge(int new_id,node node_g,node node_h,float height) {
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
    void init(const NumericMatrix X) {};
    median() : GTMethod() {};
  };
  
  class chisq : public GTMethod {
  public:
    void init(const NumericMatrix X) {
      int D = X.ncol();
      int T = sum(X);
      NumericVector wt(D);
      for(int d=0; d<D; ++d){
        wt(d)=T/sum(X(_,d));
      }
      w=wt;
    };
    float dist(node node_g,node node_h) {
      return dist_chisq(node_g.x,node_h.x,w);
    };
    node merge(int new_id,node node_g,node node_h,float height) {
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
    void init(const NumericMatrix X) {};
    node init_node(int id,NumericVector x) {
      node cnode;
      cnode.id   = id;
      cnode.x    = x; 
      cnode.size = 1;
      float ldm = log_dirichlet_multinom(x,beta);
      List cstats = List::create(Named("d",0),
                                 Named("pi",0),
                                 Named("L",ldm),
                                 Named("Lt",ldm),
                                 Named("r",0));
      cnode.stats = cstats;
      return cnode;
    };
    float dist(node node_g,node node_h) {
      float dg = node_g.stats["d"];
      float dh = node_h.stats["d"];
      float Ltg = node_g.stats["Lt"];
      float Lth = node_h.stats["Lt"];
      float dn   = ladd(lgamma(node_g.size+node_h.size),dg+dh);
      float pi   = lgamma(node_g.size+node_h.size)-dn;
      NumericVector x   = node_g.x + node_h.x;
      float L    = log_dirichlet_multinom(x,beta);
      float Lt   = ladd(pi+L,-pi+Ltg+Lth);
      // reverse for min heap
      return Lt-pi-L;
    };
    node merge(int new_id,node node_g,node node_h,float height) {
      node new_node;
      new_node.id = new_id;
      new_node.x  = node_g.x + node_h.x;
      new_node.size= node_h.size+node_g.size;
      new_node.height = height;
      
      float dg = node_g.stats["d"];
      float dh = node_h.stats["d"];
      float Ltg = node_g.stats["Lt"];
      float Lth = node_h.stats["Lt"];
      
      float d   = ladd(lgamma(node_g.size+node_h.size),dg+dh);
      float pi  = lgamma(node_g.size+node_h.size)-d;
      float L   = log_dirichlet_multinom(new_node.x,beta);
      float Lt  = ladd(pi+L,-pi+Ltg+Lth);
      float r   = pi+L-Lt;
      
      List cstats = List::create(Named("d",d),
                                 Named("pi",pi),
                                 Named("L",L),
                                 Named("Lt",Lt),
                                 Named("r",r));
      new_node.stats = cstats;
      return new_node;
    }
    bayes_mom(float beta_val) {
      beta = beta_val;
    };
  private:
    float beta;
  };

}


#endif
