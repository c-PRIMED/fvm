#ifndef _DISTFUNCTFIELDS_H_
#define _DISTFUNCTFIELDS_H_


#include "Field.h"
#include "DistFunctModel.h"

struct DistFunctFields
{
  DistFunctFields(const string baseName);

  Field velocity;
  Field pressure;
  Field massFlux;
  Field velocityGradient;
  Field pressureGradient;
  Field momentumFlux;
  Field viscosity;
  Field density;
  Field continuityResidual;
  Field velocityN1;
  Field velocityN2;
};

#endif
