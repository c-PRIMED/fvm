#include "ThermalFields.h"

ThermalFields::ThermalFields(const string baseName) :
  temperature(baseName + ".temperature"),
  heatFlux(baseName + ".heatFlux"),
  temperatureGradient(baseName + ".temperatureGradient"),
  conductivity(baseName + ".conductivity")
{}

