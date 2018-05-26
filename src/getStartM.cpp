#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getStartM(NumericVector N, NumericVector w,
                        double delta, double dt, double Ystar=3.5) {

  int n = N.length();
  NumericMatrix M(n, n);
  double Nj  = 0.0;
  double Cij = 0.0;

  for(int i=0; i < n; i++)
  {
    double Ysi = Ystar*w[i]*N[i]*dt;
    for(int j=0; j < n; j++)
    {
      Nj  = floor(N[j]);
      Cij = floor(std::min(Ysi/w[j], delta*Nj));
      if(Cij > 1) M(i,j) = log(Nj/(Nj-Cij));
    }
  }

  return M;

}



