// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <string>

using namespace std;

#include "PC.h"


typedef PC<3,2> PCType;

int main(int argc, char *argv[])
{
  PCType a;
  a = 0.1, 0.01;
  
  PCType b = a*a;
  PCType c = a+3;
  PCType d = a*4;


  cout << "mean value of a " << a._data[0] << endl;
  cout << "std dev of a " << stdDev(a) << endl;
  
  cout << "mean value of b " << b._data[0] << endl;
  cout << "std dev of b " << stdDev(b) << endl;

  cout << "mean value of c " << c._data[0] << endl;
  cout << "std dev of c " << stdDev(c) << endl;

  cout << "mean value of d " << d._data[0] << endl;
  cout << "std dev of d " << stdDev(d) << endl;

  return 0;
}
