// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GEOMFIELDS_H_
#define _GEOMFIELDS_H_

#include "Field.h"

#ifndef SWIG
class StorageSite;
class Matrix;
#endif

struct GeomFields
{
  GeomFields(const string baseName);

  Field coordinate;
  Field coordinateN1;
  Field coordinate0; // value at t=0
  Field coordinateK1; //value at previous iteration
  Field area;
  Field areaN1;
  Field areaMag;
  Field volume;
  Field volumeN1;
  Field volumeN2;
  Field sweptVolDot;
  Field sweptVolDotN1;
  Field gridFlux;
  Field faceVel;
  Field nodeDisplacement;
  Field nodeDisplacementN1;
  Field boundaryNodeNormal;
  Field dirichletNodeDisplacement;
  Field displacementOptions;
  Field ibType;
  Field ibTypeN1;
  Field ibFaceIndex;
  Field fineToCoarse;
  Field finestToCoarse;
  // this file gets directly included in a swig ineterface definition
  // file hence protect the following
#ifndef SWIG

  
  typedef pair<const StorageSite*, const StorageSite*> SSPair;

  mutable map<SSPair,shared_ptr<Matrix> > _interpolationMatrices;
  
#endif

};

#endif
