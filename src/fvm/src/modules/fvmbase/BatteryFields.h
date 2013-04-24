// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BATTERYFIELDS_H_
#define _BATTERYFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct BatterySpeciesFields
{
  BatterySpeciesFields(const string baseName);

  Field concentration;
  Field massFlux;
  Field diffusivity;
  Field convectionFlux;
  Field concentrationN1;
  Field concentrationN2;
  Field one;                      //used to fill in density
};

struct BatteryModelFields
{
  BatteryModelFields(const string baseName);

  Field potential; 
  Field potentialN1;
  Field potential_flux;  
  Field potential_gradient;
  Field conductivity;
  Field lnLithiumConcentration;
  Field lnLithiumConcentrationGradient;
  Field speciesGradient;  //  used temporarily for each species linearization
  Field temperature;
  Field temperatureN1;
  Field temperatureN2;
  Field rhoCp;
  Field heatFlux;
  Field temperatureGradient;
  Field thermalConductivity;
  Field heatSource;
  Field potentialSpeciesTemp;
  Field potentialSpeciesTempN1;
  Field potentialSpeciesTempN2;
  Field potentialSpeciesTempFlux;
  Field potentialSpeciesTempGradient;
  Field potentialSpeciesTempDiffusivity;
};

#endif
