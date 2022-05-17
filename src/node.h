#ifndef NODE
#define NODE

#include <Rcpp.h>

using namespace Rcpp;


struct node
{
  int id;
  int size;
  List stats;
  NumericVector x;
  double height;
  std::map<int,double,std::greater<double>> neibs;
};

#endif
