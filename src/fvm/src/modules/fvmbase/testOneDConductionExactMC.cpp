// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <iostream>
#include <fstream>

#include <math.h>
using namespace std;
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_histogram.h>

#define NCELLS 10
#define NTRIES 100000
int main()
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r,455465636);
  
  double sigma1=0.1;
  
  double **Tc = new double*[NCELLS];
  double *xc = new double[NCELLS];
  for(int c=0; c<NCELLS; c++)
      Tc[c] = new double[NTRIES];
  
  double dx = 1.0/NCELLS;

  xc[0] = dx/2.0;
  for(int c=1; c<NCELLS; c++)
    xc[c] = xc[c-1] + dx;

  // increment inside the loop
  for(int i=0; i<NTRIES; )
  {
      double e = (gsl_ran_ugaussian(r)*sigma1 +0.);
      if (e<-1.0) continue;
      double ln_term = log(1 + e);
      for(int c=0; c<NCELLS; c++)
      {
          if (e == 0.0)
          {
              Tc[c][i] = xc[c];
          }
          else
          {
              Tc[c][i] = log(1+e*xc[c])/ln_term;
          }
      }
      i++;
  }
  

  ofstream mean_file, sd_file;
  mean_file.open("/tmp/mean-mc.dat");
  sd_file.open("/tmp/sd-mc.dat");
  
  for(int c=0; c<NCELLS; c++)
  {
      double Tmean = gsl_stats_mean(Tc[c],1,NTRIES);
      double Tvariance = gsl_stats_variance_m(Tc[c],1,NTRIES,Tmean);
      mean_file << xc[c] << "  " << Tmean << endl;
      sd_file << xc[c] << "  " << sqrt(Tvariance) << endl;
  }
  mean_file.close();
  sd_file.close();
  exit(0);
}
