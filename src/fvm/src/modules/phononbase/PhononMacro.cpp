#include "PhononMacro.h"

PhononMacro::PhononMacro(const std::string baseName) :
  temperature(baseName + ".temperature"),
  deltaT(baseName + ".deltaT"),
  e0(baseName + ".e0"),
  TlResidual(baseName + "TlResidual"),
  TlInjected(baseName + "TlInjected"),
  TlFASCorrection(baseName + "TlFASCorrection"),
  heatFlux(baseName + ".heatFlux"),
  lam(baseName + ".lam"),
  zero(baseName + "zero"),
  one(baseName + "one")
{}
