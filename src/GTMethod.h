#ifndef GTMETHOD
#define GTMETHOD

#include <RcppArmadillo.h>
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

  
  
class bayes_dcsbm {
  public:
    void init() {
      
    };
  bayesian_dcsbm_node init_node(int id,bayesian_dcsbm_node gnode) {
      bayesian_dcsbm_node cnode;
      cnode.id   = id;
      cnode.size = 1;
      cnode.intra_e = gnode.intra_e;
      cnode.out_edges = gnode.out_edges;
      cnode.in_edges= gnode.in_edges;
      cnode.in_deg=gnode.in_deg;
      cnode.out_deg=gnode.out_deg;
      cnode.height = 0;
      return cnode;
    };
  double dist(bayesian_dcsbm_node * node_g,bayesian_dcsbm_node * node_h,std::set<int> active_nodes,std::vector<bayesian_dcsbm_node> graph) {
    
    if(node_g->id>node_h->id){
      bayesian_dcsbm_node * nt;
      nt = node_g;
      node_g=node_h;
      node_h=nt;
    }
    // calcul de l'icl 
    // before merge
    int g_in_deg = node_g->in_deg;
    int h_in_deg = node_h->in_deg;
    int g_out_deg =node_g->out_deg;
    int h_out_deg = node_h->out_deg;
    int msize = node_g->size+node_h->size;
    int m_out_deg = h_out_deg+g_out_deg;
    int m_in_deg = h_in_deg+g_in_deg;
    // degree correction
    double icl_before = lgamma(node_g->size)+ g_in_deg*log(node_g->size)-lgamma(node_g->size+g_in_deg);
    icl_before += lgamma(node_h->size)+ h_in_deg*log(node_h->size)-lgamma(node_h->size+h_in_deg);
    icl_before += lgamma(node_g->size)+ g_out_deg*log(node_g->size)-lgamma(node_g->size+g_out_deg);
    icl_before += lgamma(node_h->size)+ h_out_deg*log(node_h->size)-lgamma(node_h->size+h_out_deg);
    double icl_after =  lgamma(msize)+ m_in_deg*log(msize)-lgamma(msize+m_in_deg);
    icl_after += lgamma(msize)+ m_out_deg*log(msize)-lgamma(msize+m_out_deg);
    // besoin des tailles de tous les clusters actifs
    std::map<int, int> g_out_edges = node_g->out_edges;
    std::map<int, int> h_out_edges = node_h->out_edges;
    std::map<int, int> g_in_edges = node_g->in_edges;
    std::map<int, int> h_in_edges = node_h->in_edges;;
    for(auto it=active_nodes.begin(); it !=active_nodes.end(); ++it){
      int cs = graph[*it].size;
      double nbl_og = 0;
      double nbl_oh = 0;
      double nbl_ig = 0;
      double nbl_ih = 0;
      
      if(*it!=node_g->id){
        auto it_og =g_out_edges.find(*it);
        if ( it_og != g_out_edges.end()){
          nbl_og=it_og->second;
          icl_before+=lgamma(nbl_og+1)+nbl_og*log(lambda_out);
        }
        icl_before -= (nbl_og+1)*log(lambda_out*node_g->size*cs+1);
      }
      
      if(*it!=node_h->id){
        auto it_oh =h_out_edges.find(*it);
        if ( it_oh != h_out_edges.end()){
          nbl_oh=it_oh->second;
          icl_before+=lgamma(nbl_oh+1)+nbl_oh*log(lambda_out);
        }
        icl_before -= (nbl_oh+1)*log(lambda_out*node_h->size*cs+1);
      }
      
      if(*it!=node_h->id && *it!=node_g->id){
        
        auto it_ig =g_in_edges.find(*it);
        
        if ( it_ig != g_in_edges.end()){
          nbl_ig=it_ig->second;
          icl_before+=lgamma(nbl_ig+1)+nbl_ig*log(lambda_out);
        }
        icl_before -= (nbl_ig+1)*log(lambda_out*node_g->size*cs+1);
        
        auto it_ih =h_in_edges.find(*it);
        
        if ( it_ih != h_in_edges.end()){
          nbl_ih=it_ih->second;
          icl_before+=lgamma(nbl_ih+1)+nbl_ih*log(lambda_out);
        }
        icl_before -= (nbl_ih+1)*log(lambda_out*node_h->size*cs+1);
        
        icl_after += lgamma(nbl_oh+nbl_og+1)+(nbl_oh+nbl_og)*log(lambda_out);
        icl_after -= (nbl_oh+nbl_og+1)*log(lambda_out*msize*cs+1);          
        
        icl_after += lgamma(nbl_ih+nbl_ig+1)+(nbl_ih+nbl_ig)*log(lambda_out);
        icl_after -= (nbl_ih+nbl_ig+1)*log(lambda_out*msize*cs+1);
        
      }
      
    }
    
    
    int g_intra = node_g->intra_e;
    icl_before += lgamma(g_intra+1)+(g_intra)*log(lambda_in);
    icl_before -= (g_intra+1)*log(lambda_in*node_g->size*node_g->size+1);
    
    
    int h_intra = node_h->intra_e;
    icl_before += lgamma(h_intra+1)+(h_intra)*log(lambda_in);
    icl_before -= (h_intra+1)*log(lambda_in*node_h->size*node_h->size+1);
    
    int nbl_gh = 0;
    auto it_exist_hg =h_out_edges.find(node_g->id);
    if ( it_exist_hg != h_out_edges.end()){
      nbl_gh+=it_exist_hg->second;
    }
    auto it_exist_gh =g_out_edges.find(node_h->id);
    if ( it_exist_gh != g_out_edges.end()){
      nbl_gh+=it_exist_gh->second;
    }

    icl_after += lgamma(g_intra+h_intra+nbl_gh+1)+(g_intra+h_intra+nbl_gh)*log(lambda_in);
    icl_after -= (g_intra+h_intra+nbl_gh+1)*log(lambda_in*msize*msize+1);
    
    return icl_before-icl_after;
    
    
  };
  bayesian_dcsbm_node merge(int new_id, bayesian_dcsbm_node * node_g, bayesian_dcsbm_node * node_h,double height) {
    bayesian_dcsbm_node  new_node;
    new_node.id = new_id;
    new_node.size= node_h->size+node_g->size;

    new_node.in_deg = node_h->in_deg+node_g->in_deg;
    new_node.out_deg = node_h->out_deg+node_g->out_deg;
    
    // attention les liens g->h et h->g deviennent des lien internes  
    new_node.intra_e = node_h->intra_e+node_g->intra_e;
    std::map<int, int> h_out_edges = node_h->out_edges;
    auto it_exist_hg =h_out_edges.find(node_g->id);
    if ( it_exist_hg != h_out_edges.end()){
      new_node.intra_e=new_node.intra_e+it_exist_hg->second;
      h_out_edges.erase(node_g->id);
    }
    std::map<int, int> g_out_edges = node_g->out_edges;
    auto it_exist_gh =g_out_edges.find(node_h->id);
    if ( it_exist_gh != g_out_edges.end()){
      new_node.intra_e=new_node.intra_e+it_exist_gh->second;
      g_out_edges.erase(node_h->id);
    }
    std::map<int, int> out_e = merge_edges(g_out_edges,h_out_edges);
    
    std::map<int, int> g_in_edges = node_g->in_edges;
    g_in_edges.erase(node_h->id);
    std::map<int, int> h_in_edges = node_h->in_edges;
    h_in_edges.erase(node_g->id); 
    std::map<int, int> in_e = merge_edges(g_in_edges,h_in_edges);
    
    new_node.in_edges=in_e;
    new_node.out_edges=out_e;
    return new_node;
  }
  bayes_dcsbm(double lambda_in_val,double lambda_out_val) {
    lambda_in= lambda_in_val;
    lambda_out= lambda_out_val;
  };
  private:
    double lambda_in;
    double lambda_out;
  };
  
  
  
  
  class bayes_pdcsbm {
  public:
    void init() {
      
    };
    bayesian_dcsbm_node init_node(int id,bayesian_dcsbm_node gnode) {
      bayesian_dcsbm_node cnode;
      cnode.id   = id;
      cnode.size = 1;
      cnode.intra_e = gnode.intra_e;
      cnode.out_edges = gnode.out_edges;
      cnode.in_edges= gnode.in_edges;
      cnode.in_deg=gnode.in_deg;
      cnode.out_deg=gnode.out_deg;
      cnode.height = 0;
      return cnode;
    };
    double dist(bayesian_dcsbm_node * node_g,bayesian_dcsbm_node * node_h,int nb_externe_pairs,int nb_externe_links) {
      
      if(node_g->id>node_h->id){
        bayesian_dcsbm_node * nt;
        nt = node_g;
        node_g=node_h;
        node_h=nt;
      }
      // calcul de l'icl 
      // before merge
      int g_in_deg = node_g->in_deg;
      int h_in_deg = node_h->in_deg;
      int g_out_deg =node_g->out_deg;
      int h_out_deg = node_h->out_deg;
      int msize = node_g->size+node_h->size;
      int m_out_deg = h_out_deg+g_out_deg;
      int m_in_deg = h_in_deg+g_in_deg;
      // degree correction
      double icl_before = lgamma(node_g->size)+ g_in_deg*log(node_g->size)-lgamma(node_g->size+g_in_deg);
      icl_before += lgamma(node_h->size)+ h_in_deg*log(node_h->size)-lgamma(node_h->size+h_in_deg);
      icl_before += lgamma(node_g->size)+ g_out_deg*log(node_g->size)-lgamma(node_g->size+g_out_deg);
      icl_before += lgamma(node_h->size)+ h_out_deg*log(node_h->size)-lgamma(node_h->size+h_out_deg);
      double icl_after =  lgamma(msize)+ m_in_deg*log(msize)-lgamma(msize+m_in_deg);
      icl_after += lgamma(msize)+ m_out_deg*log(msize)-lgamma(msize+m_out_deg);
      
      
      // besoin des tailles de tous les clusters actifs
      std::map<int, int> g_out_edges = node_g->out_edges;
      std::map<int, int> h_out_edges = node_h->out_edges;
      std::map<int, int> g_in_edges = node_g->in_edges;
      std::map<int, int> h_in_edges = node_h->in_edges;;


      icl_before += lgamma(nb_externe_links+1)+(nb_externe_links)*log(lambda_ext);
      icl_before -= (nb_externe_links+1)*log(lambda_ext*nb_externe_pairs+1);
      
      
      int g_intra = node_g->intra_e;
      icl_before += lgamma(g_intra+1)+(g_intra)*log(lambda_in);
      icl_before -= (g_intra+1)*log(lambda_in*node_g->size*node_g->size+1);
      
      
      int h_intra = node_h->intra_e;
      icl_before += lgamma(h_intra+1)+(h_intra)*log(lambda_in);
      icl_before -= (h_intra+1)*log(lambda_in*node_h->size*node_h->size+1);
      
      int nbl_gh = 0;
      auto it_exist_hg =h_out_edges.find(node_g->id);
      if ( it_exist_hg != h_out_edges.end()){
        nbl_gh+=it_exist_hg->second;
      }

      auto it_exist_gh =g_out_edges.find(node_h->id);
      if ( it_exist_gh != g_out_edges.end()){
        nbl_gh+=it_exist_gh->second;
      }

      icl_after += lgamma(g_intra+h_intra+nbl_gh+1)+(g_intra+h_intra+nbl_gh)*log(lambda_in);
      icl_after -= (g_intra+h_intra+nbl_gh+1)*log(lambda_in*msize*msize+1);
      
      int new_nb_externe_links= nb_externe_links-nbl_gh;
      int new_nb_externe_pairs = nb_externe_pairs -2*node_g->size*node_h->size;
      icl_after += lgamma(new_nb_externe_links+1)+(new_nb_externe_links)*log(lambda_ext);
      icl_after -= (new_nb_externe_links+1)*log(lambda_ext*new_nb_externe_pairs+1);

      
      return icl_before-icl_after;
      
    };
    bayesian_dcsbm_node merge(int new_id, bayesian_dcsbm_node * node_g, bayesian_dcsbm_node * node_h,double height) {
      bayesian_dcsbm_node  new_node;
      new_node.id = new_id;
      new_node.size= node_h->size+node_g->size;
      
      new_node.in_deg = node_h->in_deg+node_g->in_deg;
      new_node.out_deg = node_h->out_deg+node_g->out_deg;
      
      // attention les liens g->h et h->g deviennent des lien internes  
      new_node.intra_e = node_h->intra_e+node_g->intra_e;
      std::map<int, int> h_out_edges = node_h->out_edges;
      auto it_exist_hg =h_out_edges.find(node_g->id);
      if ( it_exist_hg != h_out_edges.end()){
        new_node.intra_e=new_node.intra_e+it_exist_hg->second;
        h_out_edges.erase(node_g->id);
      }
      std::map<int, int> g_out_edges = node_g->out_edges;
      auto it_exist_gh =g_out_edges.find(node_h->id);
      if ( it_exist_gh != g_out_edges.end()){
        new_node.intra_e=new_node.intra_e+it_exist_gh->second;
        g_out_edges.erase(node_h->id);
      }
      std::map<int, int> out_e = merge_edges(g_out_edges,h_out_edges);
      
      std::map<int, int> g_in_edges = node_g->in_edges;
      g_in_edges.erase(node_h->id);
      std::map<int, int> h_in_edges = node_h->in_edges;
      h_in_edges.erase(node_g->id); 
      std::map<int, int> in_e = merge_edges(g_in_edges,h_in_edges);
      
      new_node.in_edges=in_e;
      new_node.out_edges=out_e;
      return new_node;
    }
    bayes_pdcsbm(double lambda_val_in,double lambda_val_ext) {
      lambda_in= lambda_val_in;
      lambda_ext= lambda_val_ext;
    };
  private:
    double lambda_in;
    double lambda_ext;
  };
  
}





#endif
