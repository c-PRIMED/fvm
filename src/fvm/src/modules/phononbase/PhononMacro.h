#ifndef _PHONONMACRO_H_
#define _PHONONMACRO_H_

#include "Field.h"

using namespace std;

struct PhononMacro
{

  typedef shared_ptr<Field> FieldPtr;
  typedef vector<FieldPtr> FieldVector;
  typedef map<int,FieldVector*> FieldVectorMap;

  PhononMacro(const string baseName);

  Field temperature;
  Field deltaT;
  Field e0;
  Field TlResidual;
  Field TlInjected;
  Field TlFASCorrection;
  Field heatFlux;
  Field lam;
  FieldVectorMap BranchTemperatures;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density

  Field& getModeTemp(int mesh, int mode);

};

#endif
