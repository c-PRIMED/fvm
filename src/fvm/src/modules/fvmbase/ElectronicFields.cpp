#include "ElectronicFields.h"

ElectronicFields::ElectronicFields(const string baseName) :
  potential(baseName + ".potential"),
  potentialFlux(baseName + ".potentialFlux"),
  potentialGradient(baseName + ".potentialGradient"),
  permittivity(baseName + ".permittivity"),
  source(baseName + ".source"),

  charge(baseName + "charge"),
  chargeFlux(baseName + "chargeFlux"),
  diffusivity(baseName + "diffusivity"),
  convectionFlux(baseName + "convectionFlux"),
  chargeGradient(baseName + "chargeGradient"),
  chargeN1(baseName + "chargeN1"), 
  chargeN2(baseName + "chargeN2"),

  tunnelingCharge(baseName+"tunnelingCharge"),
  
  zero(baseName + "zero"),
  one(baseName + "one")
{}

