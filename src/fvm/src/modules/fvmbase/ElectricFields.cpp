#include "ElectricFields.h"

ElectricFields::ElectricFields(const string baseName) :
  potential(baseName + ".potential"),
  potential_flux(baseName + ".potential_flux"),
  electric_field(baseName + ".electric_field"),
  dielectric_constant(baseName + ".dielectric_constant"),
  total_charge(baseName + ".total_charge"),

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

