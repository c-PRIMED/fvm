#include "FlowFields.h"

FlowFields::FlowFields(const string baseName) :
  velocity(baseName + ".velocity"),
  pressure(baseName + ".pressure"),
  massFlux(baseName + ".massFlux"),
  velocityGradient(baseName + ".velocityGradient"),
  pressureGradient(baseName + ".pressureGradient"),
  momentumFlux(baseName + ".momentumFlux"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  continuityResidual(baseName + ".continuityResidual"),
  velocityN1(baseName + ".velocityN1"),
  velocityN2(baseName + ".velocityN2"),
  stress(baseName + ".stress")
{}

