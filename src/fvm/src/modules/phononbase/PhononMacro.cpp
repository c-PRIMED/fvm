#include "PhononMacro.h"

PhononMacro::PhononMacro(const std::string baseName) :
  temperature(baseName + ".temperature"),
  heatFlux(baseName + ".heatFlux"),
  zero(baseName + "zero"),
  one(baseName + "one")
{}
