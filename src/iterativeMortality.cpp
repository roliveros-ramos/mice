#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List iterativeMortality(NumericMatrix Mj, NumericVector add, NumericVector w,
                        NumericVector N, NumericVector Yso, int niter)
{

  //NumericMatrix Mj(rMj);
  //NumericVector add(radd);
  //NumericVector w(rw);
  //NumericVector N(rN);
  //NumericVector Yso(rYso);
  //int niter = rniter;

  int nrow = Mj.nrow();
  int ncol = Mj.ncol();

  NumericMatrix Cj(nrow, ncol);
  NumericVector Yj(ncol);

    for (int k=0; k<niter; k++)
    {
      NumericVector Zj(nrow);
      NumericVector Dj(nrow);
      for(int i=0; i<nrow; i++)
      {
        for(int j=0; j<ncol; j++)
        {
          Zj[i] += Mj(i,j);
        }
        Zj[i] += add[i];
        Dj[i] = 1 - exp(-Zj[i]);
      }

      NumericVector cfj(ncol, 1.0);
      double ans = 0.0;
      for(int i=0; i<nrow; i++)
      {
        for(int j=0; j<ncol; j++)
        {
          ans = (Mj(i,j)/Zj[i])*Dj[i]*N[i];
          if(traits::is_nan<REALSXP>(ans)) {
            Cj(i,j) = 0;
          } else {
            Cj(i,j) = ans;
          }
          Yj[j] += w[i]*Cj(i,j);
        }
      }

      double cfj_tmp = 0.0;
      for(int j=0; j<ncol; j++)
      {
        cfj_tmp = Yso[j]/Yj[j];
        if(traits::is_nan<REALSXP>(cfj_tmp)) cfj_tmp = 0;
        if(cfj_tmp < 1) cfj[j] = cfj_tmp;
      }

      for(int i=0; i<nrow; i++)
      {
        for(int j=0; j<ncol; j++)
        {
          Mj(i,j) = cfj[j]*Mj(i,j);
        }
      }

    }

    return List::create(Named("M", Mj), Named("Yj", Yj));

}

