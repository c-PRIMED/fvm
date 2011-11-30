#ifndef _PHONONMACRO_H_
#define _PHONONMACRO_H_

#include "Field.h"

using namespace std;

struct PhononMacro
{
  PhononMacro(const string baseName);

  Field temperature;
  Field deltaT;
  Field e0;
  Field e0Residual;
  Field e0Injected;
  Field e0FASCorrection;
  Field heatFlux;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
};

#endif
