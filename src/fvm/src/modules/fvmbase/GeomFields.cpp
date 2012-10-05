// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "GeomFields.h"

GeomFields::GeomFields(const string baseName) :
  coordinate(baseName+"coordinate"),
  coordinateN1(baseName+"coordinateN1"),
  coordinate0(baseName+"coordinate0"),
  coordinateK1(baseName+"coordinateK1"),
  area(baseName+"area"),
  areaN1(baseName+"areaN1"),
  areaMag(baseName+"areaMag"),
  volume(baseName+"volume"),
  volumeN1(baseName+"volumeN1"),
  volumeN2(baseName+"volumeN2"),
  sweptVolDot(baseName+"sweptVolDot"),
  sweptVolDotN1(baseName+"sweptVolDotN1"),
  gridFlux(baseName+"gridFlux"),
  faceVel(baseName+"faceVel"),
  nodeDisplacement(baseName+"nodeDisplacement"),
  nodeDisplacementN1(baseName+"nodeDisplacementN1"),
  boundaryNodeNormal(baseName+"boundaryNodeNormal"),
  dirichletNodeDisplacement(baseName+"dirichletNodeDisplacement"),
  displacementOptions(baseName+"displacementOptions"),
  ibType(baseName+"ibType"),
  ibTypeN1(baseName+"ibTypeN1"),
  ibFaceIndex(baseName+"ibFaceIndex"),
  fineToCoarse(baseName+"fineToCoarse"),
  finestToCoarse(baseName+"finestToCoarse")
{}

