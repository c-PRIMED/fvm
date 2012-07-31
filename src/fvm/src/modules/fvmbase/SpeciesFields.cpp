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
  elecPotential(baseName + "elecPotential")
{}

SpeciesModelFields::SpeciesModelFields(const string baseName) :
  speciesGradient(baseName + ".speciesGradient")
{}
