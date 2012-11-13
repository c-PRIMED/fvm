// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "BatteryFields.h"

BatterySpeciesFields::BatterySpeciesFields(const string baseName) :
  massFraction(baseName + ".massFraction"),
  massFlux(baseName + ".massFlux"),
  diffusivity(baseName + ".diffusivity"),
  convectionFlux(baseName + ".convectionFlux"),
  massFractionN1(baseName + ".massFractionN1"),
  massFractionN2(baseName + ".massFractionN2"),
  one(baseName + "one")
{}

BatteryModelFields::BatteryModelFields(const string baseName) :
  potential(baseName + "potential"),
  potential_flux(baseName + ".potential_flux"),
  potential_gradient(baseName + ".potential_gradient"),
  conductivity(baseName + ".conductivity"),
  lnLithiumConcentration(baseName + "lnLithiumConcentration"),
  lnLithiumConcentrationGradient(baseName + "lnLithiumConcentrationGradient"),
  speciesGradient(baseName + ".speciesGradient"),
  potentialAndSpecies(baseName + ".potentialAndSpecies"),
  potentialAndSpeciesN1(baseName + ".potentialAndSpeciesN1"),
  potentialAndSpeciesN2(baseName + ".potentialAndSpeciesN2"),
  potentialAndSpeciesFlux(baseName + ".potentialAndSpeciesFlux"),
  potentialAndSpeciesGradient(baseName + ".potentialAndSpeciesGradient"),
  potentialAndSpeciesDiffusivity(baseName + ".potentialAndSpeciesDiffusivity")
{}
