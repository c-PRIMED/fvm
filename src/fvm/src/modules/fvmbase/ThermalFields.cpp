#include "ThermalFields.h"

ThermalFields::ThermalFields(const string baseName) :
  temperature(baseName + ".temperature"),
  temperatureN1(baseName + ".temperatureN1"),
  temperatureN2(baseName + ".temperatureN2"),
  heatFlux(baseName + ".heatFlux"),
  temperatureGradient(baseName + ".temperatureGradient"),
  conductivity(baseName + ".conductivity"),
  source(baseName + ".source"),
  convectionFlux(baseName + ".convectionFlux"),
  zero(baseName + "zero"),
  one(baseName + "one"),
  specificHeat(baseName + ".specificHeat")
{}


