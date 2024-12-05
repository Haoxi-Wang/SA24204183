#include <Rcpp.h>
using namespace Rcpp;
//' @title Rcpp test
//' @description test
//' @param x int
//' @param y int
//' @return max \code{n}
//' @export
// [[Rcpp::export]]
int max(int x, int y) {
  if (x > y) {
    return x;
  } else if (x < y) {
    return y;  
  } else {
    return x;  
  }
}
  