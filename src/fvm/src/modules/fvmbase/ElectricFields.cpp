// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "ElectricFields.h"

ElectricFields::ElectricFields(const string baseName) :
  potential(baseName + ".potential"),
  potential_flux(baseName + ".potential_flux"),
  potential_gradient(baseName + ".potential_gradient"),
  electric_field(baseName + ".electric_field"),
  dielectric_constant(baseName + ".dielectric_constant"),
  total_charge(baseName + ".total_charge"),
  conduction_band(baseName + ".conduction_band"),
  valence_band(baseName + ".valence_band"),
  electron_totaltraps(baseName + ".electron_totaltraps"),
  free_electron_capture_cross(baseName + ".free_electron_capture_cross"), 
  transmission(baseName + ".transmission"),
  electron_velocity(baseName + ".electron_velocity"),
  charge(baseName + "charge"),
  chargeFlux(baseName + "chargeFlux"),
  diffusivity(baseName + "diffusivity"),
  convectionFlux(baseName + "convectionFlux"),
  chargeGradient(baseName + "chargeGradient"),
  chargeN1(baseName + "chargeN1"), 
  chargeN2(baseName + "chargeN2"),
  zero(baseName + "zero"),
  one(baseName + "one"),
  //oneD_column(baseName + "oneD_column"),
  force(baseName + "force"),
  speciesConcentration(baseName + "speciesConcentration"),
  lnSpeciesConcentration(baseName + "lnSpeciesConcentration"),
  lnSpeciesConcentrationGradient(baseName + "lnSpeciesconCentrationGradient")
  
{}

