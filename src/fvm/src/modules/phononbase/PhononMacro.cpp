#include "PhononMacro.h"

PhononMacro::PhononMacro(const std::string baseName) :
  temperature(baseName + ".temperature"),
  deltaT(baseName + ".deltaT"),
  e0(baseName + ".e0"),
  e0Residual(baseName + "e0Residual"),
  e0Injected(baseName + "e0Injected"),
  e0FASCorrection(baseName + "e0FASCorrection"),
  heatFlux(baseName + ".heatFlux"),
  lam(baseName + ".lam"),
  zero(baseName + "zero"),
  one(baseName + "one")
{}
