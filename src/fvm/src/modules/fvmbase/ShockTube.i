#ifndef _SHOCKTUBE_H_
#define _SHOCKTUBE_H_

%{
#include "ShockTube.h"
%}

using namespace std;

template<class T>
class ShockTube
{
public:
  ShockTube(const int nCells);

  ArrayBase& getSolution();
  ArrayBase& getPressure();
  ArrayBase& getDensity();
  ArrayBase& getVelocity();
  ArrayBase& getMachNumber();
  
  void solve(const double dt, const int nsteps);
  void init();
  T pL;
  T pR;
  T rhoL;
  T rhoR;
  T uL;
  T uR;
};


%template(ShockTubeA) ShockTube< ATYPE_STR >;

#endif
