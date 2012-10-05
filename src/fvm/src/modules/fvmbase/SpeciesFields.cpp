// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "SpeciesFields.h"

SpeciesFields::SpeciesFields(const string baseName) :
  massFraction(baseName + ".massFraction"),
  massFlux(baseName + ".massFlux"),
  diffusivity(baseName + ".diffusivity"),
  source(baseName + ".source"),
  convectionFlux(baseName + ".convectionFlux"),
  massFractionN1(baseName + ".massFractionN1"),
  massFractionN2(baseName + ".massFractionN2"),
  zero(baseName + "zero"),
  one(baseName + "one"),
  // For coupling to ElectricModel
  elecPotential(baseName + "elecPotential"),
  massFractionElectricModel(baseName + ".massFractionElectricModel") // mass Fraction that the Electric model last used
{}

SpeciesModelFields::SpeciesModelFields(const string baseName) :
  speciesGradient(baseName + ".speciesGradient")
{}
