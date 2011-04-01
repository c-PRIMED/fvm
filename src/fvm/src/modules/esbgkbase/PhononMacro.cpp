#include "PhononMacro.h"

PhononMacro::PhononMacro(const std::string baseName) :
  temperature(baseName + ".temperature"),
  e0(baseName + ".e0"),
  heatFlux(baseName + ".heatFlux"),
  zero(baseName + "zero"),
  one(baseName + "one")
{}
