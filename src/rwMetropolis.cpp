#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @return a random sample of size
//' @useDynLib SC19077
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100,10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]

NumericVector rwMetropolisC(double x0, int N, double sigma){
  NumericVector xC(N);
  xC[0] = x0;
  NumericVector u = runif(N);
  int k = 0;
  for(int i = 1; i < N; i++) {
    double y = (rnorm(1,xC[i-1],sigma))[0];
    double res = exp(abs(xC[i-1]) - abs(y));
    if(u[i] < res){
      xC[i] = y;
    }
    else{
      xC[i] = xC[i-1];
      k ++;}
  };
  return xC;}
  
int rw_k_C(double x0, int N, double sigma){
    NumericVector xC(N);
    xC[0] = x0;
    NumericVector u = runif(N);
    int k = 0;
    for(int i = 1; i < N; i++) {
      double y = (rnorm(1,xC[i-1],sigma))[0];
      double res = exp(abs(xC[i-1]) - abs(y));
      if(u[i] < res){
        xC[i] = y;
      }
      else{
        xC[i] = xC[i-1];
        k ++;}
    };
    return k;}
    