#ifndef _THERMALFIELDS_H_
#define _THERMALFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct ThermalFields
{
  ThermalFields(const string baseName);

  Field temperature;
  Field heatFlux;
  Field temperatureGradient;
  Field conductivity;

};

#endif
