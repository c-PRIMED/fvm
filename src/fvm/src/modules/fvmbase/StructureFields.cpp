// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "StructureFields.h"

StructureFields::StructureFields(const string baseName) :
  deformation(baseName + ".deformation"),
  elasticDeformation(baseName + ".elasticDeformation"),
  deformationGradient(baseName + ".deformationGradient"),
  deformationFlux(baseName + ".deformationFlux"),
  eta(baseName + ".eta"),
  eta1(baseName + ".eta1"),
  alpha(baseName + ".alpha"),
  density(baseName + ".density"),
  deformationN1(baseName + ".deformationN1"),
  deformationN2(baseName + ".deformationN2"),
  deformationN3(baseName + ".deformationN3"),
  tractionX(baseName + ".tractionX"),
  tractionY(baseName + ".tractionY"),
  tractionZ(baseName + ".tractionZ"),
  strainX(baseName + ".strainX"),
  strainY(baseName + ".strainY"),
  strainZ(baseName + ".strainZ"),
  plasticDiagStrain(baseName + ".plasticDiagStrain"),
  devStress(baseName + ".devStress"),
  VMStress(baseName + ".VMStress"),
  plasticStrain(baseName + ".plasticStrain"),
  temperature(baseName + ".temperature"),
  bodyForce(baseName + ".bodyForce"),
  volume0(baseName + ".volume0"),
  creepConstant(baseName + ".creepConstant")
{}

