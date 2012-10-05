// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "PlateFields.h"

PlateFields::PlateFields(const string baseName) :
  deformation(baseName + ".deformation"),
  moment(baseName + ".moment"),
  stress(baseName + ".stress"),
  residualStress(baseName + ".residualStress"),
  devStress(baseName + ".devStress"),
  VMStress(baseName + ".VMStress"),
  VMStressOut(baseName + ".VMStressOut"),
  strain(baseName + ".strain"),
  plasticStrain(baseName + ".plasticStrain"),
  plasticStrainOut(baseName + ".plasticStrainOut"),
  plasticStrainN1(baseName + ".plasticStrainN1"),
  plasticMoment(baseName + ".plasticMoment"),
  deformationGradient(baseName + ".deformationGradient"),
  deformationFlux(baseName + ".deformationFlux"),
  ym(baseName + ".ym"),
  nu(baseName + ".nu"),
  density(baseName + ".density"),
  deformationN1(baseName + ".deformationN1"),
  deformationN2(baseName + ".deformationN2"),
  deformationN3(baseName + ".deformationN3"),
  tractionX(baseName + ".tractionX"),
  tractionY(baseName + ".tractionY"),
  tractionZ(baseName + ".tractionZ"),
  thickness(baseName + ".thickness"),
  force(baseName + ".force"),
  acceleration(baseName + ".acceleration"),
  velocity(baseName + ".velocity"),
  bodyForce(baseName + ".bodyForce"),
  volume0(baseName + ".volume0")
{}

