#include <Rcpp.h>
//#include <getStartM.h>
using namespace Rcpp;

// [[Rcpp::export]]
List iterativeMortalityF(NumericVector N, NumericVector w, NumericVector delta,
                         double dt, NumericVector Ystar, NumericMatrix access,
                         NumericVector F, NumericVector add, NumericVector Mstarv,
                         double eta_crit, int niter)
{

  // initial M guess
  int n = N.length();
  NumericMatrix Mj(n, n);
  NumericVector Yso(n);
  double Ni  = 0.0;
  double dNi = 0.0;
  double Cij = 0.0;

  for(int j=0; j < n; j++)
  {
    Yso[j] = Ystar[j]*w[j]*N[j]*dt;
  }

  for(int i=0; i < n; i++)
  {
    Ni  = floor(N[i]);
    dNi = delta[i]*Ni;
    for(int j=0; j < n; j++)
    {
      Cij = floor(std::min(Yso[j]/w[i], dNi));
      if(Cij > 1) Mj(i,j) = access(i,j)*log(Ni/(Ni-Cij));
    }
  }

  // iterative solution

  NumericMatrix C(n, n);
  NumericMatrix Y(n, n);
  NumericVector Yj(n);
  NumericVector cfj(n, 1.0);
  NumericVector ratio(n, 1.0);
  NumericVector Ms(n);
  NumericVector deaths(n);

    for (int k=0; k<niter; k++)
    {
      NumericVector Z(n);

      // update of natural mortality matrix
      for(int i=0; i<n; i++)
      {
        for(int j=0; j<n; j++)
        {
          Mj(i,j) = cfj[j]*Mj(i,j);
        }
      }

      // estimation of Z and total number of deaths
      for(int i=0; i<n; i++)
      {
        for(int j=0; j<n; j++)
        {
          Z[i] += Mj(i,j);
        }

        Ms[i] = Mstarv[i]*ratio[i]; // starvation mortality
        Z[i] += (add[i] + F[i] + Ms[i])*dt;
        deaths[i] = (1 - exp(-Z[i]))*N[i]; // deaths

      }

      // estimation of predated biomass
      double ans = 0.0;
      for(int i=0; i<n; i++)
      {
        for(int j=0; j<n; j++)
        {
          ans = (Mj(i,j)/Z[i])*deaths[i];
          if(traits::is_nan<REALSXP>(ans)) {
            C(i,j) = 0;
          } else {
            C(i,j) = ans;
          }
          Y(i,j) = w[i]*C(i,j); // predated mass by prey
          Yj[j] += Y(i,j); // predated mass
        }
      }

      // updated of correction factor and starvation ratio

      double cfj_tmp = 0.0;
      double tmp = 0.0;

      for(int j=0; j<n; j++)
      {
        cfj_tmp = Yso[j]/Yj[j];
        if(traits::is_nan<REALSXP>(cfj_tmp)) cfj_tmp = 0;
        if(cfj_tmp < 1) cfj[j] = cfj_tmp;

        tmp = 1 - (Yj[j]/Yso[j])/eta_crit;
        if(traits::is_nan<REALSXP>(tmp)) tmp = 0;
        if(tmp >= 0) ratio[j] = tmp;
        if(tmp < 0) ratio[j] = 0;

      }


    } // end of k loop


    NumericVector M(n, 0.0);
    NumericVector Z(n, 0.0);
    // return M in 1/year
    for(int i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
      {
        Mj(i,j) = Mj(i,j)/dt;
        M[i] += Mj(i,j); // total natural mortality
      }
        Z[i] = M[i] + add[i] + F[i] + Ms[i]; // total mortality
    }


    return List::create(Named("M", M),
                        Named("Ms", Ms),
                        Named("Z", Z),
                        Named("F", F),
                        Named("add", add),
                        Named("preyed", Y),
                        Named("deaths", deaths),
                        Named("Mj", Mj),
                        Named("starv", ratio),
                        Named("Yj", Yj));

}

