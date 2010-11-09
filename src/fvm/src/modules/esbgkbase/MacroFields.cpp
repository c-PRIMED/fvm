#include "MacroFields.h"

MacroFields::MacroFields(const string baseName) :
  velocity(baseName + ".velocity"),
  pressure(baseName + ".pressure"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  temperature(baseName +".temperature"),
  Txx(baseName +".Txx"),
  Tyy(baseName +".Tyy"), 
  Tzz(baseName +".Tzz"), 
  Txy(baseName +".Txy"), 
  Txz(baseName +".Txz"),
  Tyz(baseName +".Tyz"),
  coeff(baseName +".coeff"),
  coeffg(baseName +".coeffg"),
  collisionFrequency(baseName+".collisionFrequency")
{}

