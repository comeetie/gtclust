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
  float height;
  std::map<int,float,std::greater<float>> neibs;
};

#endif
