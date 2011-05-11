#include "MacroFields.h"

MacroFields::MacroFields(const string baseName) :
  velocity(baseName + ".velocity"),
  pressure(baseName + ".pressure"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  temperature(baseName +".temperature"),
  collisionFrequency(baseName+".collisionFrequency"),
  Txx(baseName +".Txx"),
  Tyy(baseName +".Tyy"), 
  Tzz(baseName +".Tzz"), 
  Txy(baseName +".Txy"), 
  Txz(baseName +".Txz"),
  Tyz(baseName +".Tyz"),
  coeff(baseName +".coeff"),
  coeffg(baseName +".coeffg"),
  Entropy(baseName +"Entropy"),
  EntropyGenRate(baseName +"EntropyGenRate"), 
  EntropyGenRate_Collisional(baseName +"EntropyGenRate_Collisional"),
  force(baseName +"Force"),
  Stress(baseName +"Stress")
{}

