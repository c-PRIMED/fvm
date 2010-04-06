#include "MacroFields.h"

MacroFields::MacroFields(const string baseName) :
  velocity(baseName + ".velocity"),
  pressure(baseName + ".pressure"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  temperature(baseName +".temperature"),
  collisionFrequency(baseName+".collisionFrequency")
{}

