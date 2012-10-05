// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "FlowFields.h"

FlowFields::FlowFields(const string baseName) :
  velocity(baseName + ".velocity"),
  pressure(baseName + ".pressure"),
  massFlux(baseName + ".massFlux"),
  velocityGradient(baseName + ".velocityGradient"),
  pressureGradient(baseName + ".pressureGradient"),
  momentumFlux(baseName + ".momentumFlux"),
  viscosity(baseName + ".viscosity"),
  density(baseName + ".density"),
  continuityResidual(baseName + ".continuityResidual"),
  velocityN1(baseName + ".velocityN1"),
  velocityN2(baseName + ".velocityN2"),
  tractionX(baseName + ".tractionX"),
  tractionY(baseName + ".tractionY"),
  tractionZ(baseName + ".tractionZ"),
  stress(baseName + ".stress"),
  force(baseName + ".force"),
  eddyviscosity(baseName + ".eddyviscosity"),
  totalviscosity(baseName + ".totalviscosity"),
  utau(baseName + ".utau"),
  uparallel(baseName + ".uparallel"),
  tau(baseName + ".tau"),
  tauwall(baseName + ".tauwall"),

  //ESInterface
  InterfaceVelocity(baseName+".InterfaceVelocity"),
  InterfacePressure(baseName+".InterfacePressure"), 
  InterfaceStress(baseName+".InterfaceStress"),
  InterfaceDensity(baseName+".InterfaceDensity")

{}

