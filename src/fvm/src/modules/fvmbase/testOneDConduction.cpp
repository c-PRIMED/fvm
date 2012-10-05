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

#define NCELLS 20
int main(int argc, char *argv[])
{

  PCType kConst;
  kConst = 1.,0.1;

  OneDConduction<PCType> model(NCELLS,kConst);

  model.solve();

  Array<PCType>& Tc = *model.getSolution();

  ofstream mean_file, sd_file;
  mean_file.open("/tmp/mean-uqtk.dat");
  sd_file.open("/tmp/sd-uqtk.dat");

  double dx=1.0/NCELLS;
  double x = dx/2.0;
  for(int c=0; c<NCELLS; c++)
  {
      mean_file << x << "  " << Tc[c]._data[0] << endl;
      sd_file << x << "  " << stdDev(Tc[c]) << endl;
      x += dx;
  }

  mean_file.close();
  sd_file.close();

  return 0;
}
