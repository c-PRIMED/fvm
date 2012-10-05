// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#include "PC.h"
#include "OneDConduction.h"

typedef PC<3,1> PCType;

#define NCELLS 10
int main(int argc, char *argv[])
{

  double dkConst = 1.;
  PCType pCkConst(dkConst);
  pCkConst._data[1] = 0.1;
  
  //OneDConduction<double> dp(NCELLS,dkConst);
  OneDConduction<PCType> pCp(NCELLS,pCkConst);

  //dp.solve();
  pCp.solve();

  //Array<double>& dt = dp.getSolution();
  Array<PCType>& pCt = pCp.getSolution();

  //cout << "double solution" << endl;
  //for(int i=0; i<NCELLS; i++)
  //  cout << i << " " << dt[i] << endl;
  ofstream mean_file, variance_file;
  mean_file.open("/tmp/mean-mc.dat");
  variance_file.open("/tmp/variance-mc.dat");

  double dx=1.0/NCELLS;
  double x = dx/2.0;
  for(int c=0; c<NCELLS; c++)
  {
      double Tmean = pCt[i]._data[0] << "  "
      double T = gsl_stats_variance_m(Tc[c],1,NTRIES,Tmean);
      mean_file << x << "  " << Tmean << endl;
      variance_file << xc[c] << "  " << stdDev(pCt[i]) << endl;
  }
  
  return 0;
}
