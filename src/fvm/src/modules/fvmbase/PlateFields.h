#ifndef _PLATEFIELDS_H_
#define _PLATEFIELDS_H_


#include "Field.h"

struct PlateFields
{
  PlateFields(const string baseName);

  Field deformation;
  Field moment;
  Field deformationGradient;
  Field deformationFlux;
  Field ym;
  Field nu;
  Field density;
  Field deformationN1;
  Field deformationN2;
  Field deformationN3;
  Field tractionX;
  Field tractionY;
  Field tractionZ;
  Field thickness;
  Field force;
  Field acceleration;
  Field bodyForce;
  Field volume0;
};

#endif
