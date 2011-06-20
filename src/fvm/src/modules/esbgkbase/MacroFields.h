#ifndef _MACROPARAMETERS_H_
#define _MACROPARAMETERS_H_


#include "Field.h"

struct MacroFields
{
  MacroFields(const string baseName);
  Field velocity;
  Field pressure;
  Field viscosity;
  Field density;
  Field temperature;
  Field collisionFrequency;
  Field Txx;
  Field Tyy;
  Field Tzz;
  Field Txy;
  Field Tyz;
  Field Tzx;
  Field coeff;
  Field coeffg;
  Field Entropy;
  Field EntropyGenRate;
  Field EntropyGenRate_Collisional;
  Field force;
  Field Stress;
};

#endif
