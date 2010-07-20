#ifndef _FLOWFIELDS_H_
#define _FLOWFIELDS_H_


#include "Field.h"

struct FlowFields
{
  FlowFields(const string baseName);

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
  Field tractionX;
  Field tractionY;
  Field tractionZ;
  Field stress;
  Field force;
};

#endif
