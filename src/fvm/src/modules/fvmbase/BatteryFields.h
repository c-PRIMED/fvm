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

  Field massFraction;
  Field massFlux;
  Field diffusivity;
  Field convectionFlux;
  Field massFractionN1;
  Field massFractionN2;
  Field one;                      //used to fill in density
};

struct BatteryModelFields
{
  BatteryModelFields(const string baseName);

  Field potential;           
  Field potential_flux;  
  Field potential_gradient;
  Field conductivity;
  Field lnLithiumConcentration;
  Field lnLithiumConcentrationGradient;
  Field speciesGradient;  //  used temporarily for each species linearization
  Field potentialAndSpecies;
  Field potentialAndSpeciesN1;
  Field potentialAndSpeciesN2;
  Field potentialAndSpeciesFlux;
  Field potentialAndSpeciesGradient;
  Field potentialAndSpeciesDiffusivity;
};

#endif
