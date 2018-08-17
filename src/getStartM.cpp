#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getStartM(NumericVector N, NumericVector w,
                        NumericVector delta, double dt,
                        NumericVector Ystar, NumericMatrix access) {

  int n = N.length();
  NumericMatrix M(n, n);
  NumericVector Ys(n);
  double Ni  = 0.0;
  double dNi = 0.0;
  double Cij = 0.0;

  for(int j=0; j < n; j++)
  {
    Ys[j] = Ystar[j]*w[j]*N[j]*dt;
  }

  for(int i=0; i < n; i++)
  {
    Ni  = floor(N[i]);
    dNi = delta[i]*Ni;
    for(int j=0; j < n; j++)
    {
      Cij = floor(std::min(Ys[j]/w[i], dNi));
      if(Cij > 1) M(i,j) = access(i,j)*log(Ni/(Ni-Cij));
    }
  }

  return M;

}



