#ifndef NODE
#define NODE
#include "SuiteSparse_config/SuiteSparse_config.h"
#include "CHOLMOD/Include/cholmod.h"
#include <RcppArmadillo.h>

using namespace Rcpp;


struct abstract_node{
  int id;
  int size;
  List stats;
  NumericVector x;
  double height;
};

struct node : abstract_node
{
  std::map<int,double,std::greater<double>> neibs;
  node(){
  }
  node(abstract_node n){
    id=n.id;
    size=n.size;
    stats=n.stats;
    x=n.x;
    height=n.height;
  }
};

struct multiedge
{
  double height;
  double intra_lnbtree;
  int size;
  std::vector<std::pair<int, int>> edges;
  void add_edge(std::pair<int, int> uv){
    edges.push_back(uv);
  }
  void merge_edges(multiedge tomerge){
    size=size+tomerge.size;
    edges.insert(edges.end(),tomerge.edges.begin(),tomerge.edges.end());
  }
  void update_height(double h){
    height = h;
  }
  multiedge(int s, double h){ size=s;height=h;};
  
};


struct subfactor
{
  int * cols;
  int * nz;
  int * p;
  int * i;
  double * x;
  int size;
};


struct bayesian_node : abstract_node
{
  std::map<int,multiedge,std::greater<double>> neibs;
  std::set<int,std::greater<int>> intra_nodes;
  std::set<int> intra_pivot_edges;
  int intra_pivot;
  int i_inter;
  std::map<int, int> inter_edges;
  double lognbtree;
  bayesian_node(){
  }
  bayesian_node(abstract_node n){
    id=n.id;
    size=n.size;
    stats=n.stats;
    x=n.x;
    height=n.height;
  }
};

struct bayesian_dcsbm_node : bayesian_node{
  int intra_e=0;
  std::map<int, int> out_edges;
  std::map<int, int> in_edges;
  int in_deg=0;
  int out_deg=0;
};




#endif
