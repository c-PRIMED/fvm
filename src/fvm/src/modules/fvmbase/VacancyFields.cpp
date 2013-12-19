// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "VacancyFields.h"

VacancyFields::VacancyFields(const string baseName) :
  concentration(baseName + ".concentration"),
  concentrationN1(baseName + ".concentrationN1"),
  concentrationN2(baseName + ".concentrationN2"),
  plasticStrain(baseName + "plasticStrain"),
  vacaFlux(baseName + ".vacaFlux"),
  concentrationGradient(baseName + ".concentrationGradient"),
  concentrationGradientVector(baseName + ".concentrationGradientVector"),
  diffusioncoefficient(baseName + ".diffusioncoefficient"),
  source(baseName + ".source"),
  convectionFlux(baseName + ".convectionFlux"),
  zero(baseName + "zero"),
  one(baseName + "one"),
  specificVaca(baseName + ".specificVaca")
{}


