// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "KeFields.h"

KeFields::KeFields(const string baseName) :
  energy(baseName + ".energy"),
  energyGradient(baseName + ".energyGradient"),
  dissipation(baseName + ".dissipation"),
  dissipationGradient(baseName + ".dissipationGradient"),
  c1(baseName + ".c1"),
  c2(baseName + ".c2"),
  sigmak(baseName + ".sigmak"),
  sigmae(baseName + ".sigmae"),
  cmu(baseName + ".cmu"),
 // convectionFlux(baseName + ".convectionFlux"),
  kFlux(baseName + ".kFlux"),
  eFlux(baseName + ".eFlux"),
  energyN1(baseName + ".energyN1"),
  energyN2(baseName + ".energyN2"),
//  density(baseName + ".density"),
  dissipationN1(baseName + ".dissipationN1"),
  dissipationN2(baseName + ".dissipationN2"),
  sourcek(baseName + ".sourcek"),
  sourced(baseName + ".sourced"),
  sourcec(baseName + ".sourcec"),
  sourcep(baseName + ".sourcep")



{}



