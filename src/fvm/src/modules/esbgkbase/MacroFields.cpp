#include "MacroFields.h"

MacroFields::MacroFields(const string baseName) :
  velocity(baseName + ".velocity"),
  velocityResidual(baseName + ".velocityResidual"),
  velocityInjected(baseName + ".velocityInjected"),
  velocityFASCorrection(baseName + ".velocityFASCorrection"),
  pressure(baseName + ".pressure"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  temperature(baseName +".temperature"),
  collisionFrequency(baseName+".collisionFrequency"),
  Txx(baseName +".Txx"),
  Tyy(baseName +".Tyy"), 
  Tzz(baseName +".Tzz"), 
  Txy(baseName +".Txy"), 
  Tyz(baseName +".Tyz"),
  Tzx(baseName +".Tzx"),
  coeff(baseName +".coeff"),
  coeffg(baseName +".coeffg"),
  Entropy(baseName +"Entropy"),
  EntropyGenRate(baseName +"EntropyGenRate"), 
  EntropyGenRate_Collisional(baseName +"EntropyGenRate_Collisional"),
  force(baseName +"Force"),
  Stress(baseName +"Stress"),
  Knq(baseName+"Knq")//M300+M120+M102 for couette with uy

{}

