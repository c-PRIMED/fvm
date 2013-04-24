// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


#include "BatteryFields.h"

BatterySpeciesFields::BatterySpeciesFields(const string baseName) :
  concentration(baseName + ".concentration"),
  massFlux(baseName + ".massFlux"),
  diffusivity(baseName + ".diffusivity"),
  convectionFlux(baseName + ".convectionFlux"),
  concentrationN1(baseName + ".concentrationN1"),
  concentrationN2(baseName + ".concentrationN2"),
  one(baseName + "one")
{}

BatteryModelFields::BatteryModelFields(const string baseName) :
  potential(baseName + "potential"),
  potentialN1(baseName + "potentialN1"),
  potential_flux(baseName + ".potential_flux"),
  potential_gradient(baseName + ".potential_gradient"),
  conductivity(baseName + ".conductivity"),
  lnLithiumConcentration(baseName + "lnLithiumConcentration"),
  lnLithiumConcentrationGradient(baseName + "lnLithiumConcentrationGradient"),
  speciesGradient(baseName + ".speciesGradient"),
  temperature(baseName + ".temperature"),
  temperatureN1(baseName + ".temperatureN1"),
  temperatureN2(baseName + ".temperatureN2"),
  heatFlux(baseName + ".heatFlux"),
  temperatureGradient(baseName + ".temperatureGradient"),
  thermalConductivity(baseName + ".thermalConductivity"),
  heatSource(baseName + ".heatSource"),
  rhoCp(baseName + ".rhoCp"),
  potentialSpeciesTemp(baseName + ".potentialSpeciesTemp"),
  potentialSpeciesTempN1(baseName + ".potentialSpeciesTempN1"),
  potentialSpeciesTempN2(baseName + ".potentialSpeciesTempN2"),
  potentialSpeciesTempFlux(baseName + ".potentialSpeciesTempFlux"),
  potentialSpeciesTempGradient(baseName + ".potentialSpeciesTempGradient"),
  potentialSpeciesTempDiffusivity(baseName + ".potentialSpeciesTempDiffusivity")
{}
