// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "VacancyFields.h"

VacancyFields::VacancyFields(const string baseName) :
  concentration(baseName + ".concentration"),
  concentrationN1(baseName + ".concentrationN1"),
  concentrationN2(baseName + ".concentrationN2"),
  concentrationGradient(baseName + ".concentrationGradient"),
  diffusioncoefficient(baseName + ".diffusioncoefficient"),
  vacaFlux(baseName +".vacaFlux"),
  one(baseName + "one")
  {}


