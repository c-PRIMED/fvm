// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FractureFields.h"

FractureFields::FractureFields(const string baseName) :
  phasefieldvalue(baseName + ".phasefieldvalue"),
  phasefieldvalueN1(baseName + ".phasefieldvalueN1"),
  phasefieldvalueN2(baseName + ".phasefieldvalueN2"),
  phasefieldFlux(baseName + ".phasefieldFlux"),
  phasefieldGradient(baseName + ".phasefieldGradient"),
  conductivity(baseName + ".conductivity"),
  source(baseName + ".source"),
  sourcecoef(baseName + ".sourcecoef"),
  //convectionFlux(baseName + ".convectionFlux"),
  zero(baseName + "zero"),
  one(baseName + "one")
  //specificHeat(baseName + ".specificHeat")
{}


