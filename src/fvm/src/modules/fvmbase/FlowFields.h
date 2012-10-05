// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOWFIELDS_H_
#define _FLOWFIELDS_H_


#include "Field.h"

struct FlowFields
{
  FlowFields(const string baseName);

  Field velocity;
  Field pressure;
  Field massFlux;
  Field velocityGradient;
  Field pressureGradient;
  Field momentumFlux;
  Field viscosity;
  Field density;
  Field continuityResidual;
  Field velocityN1;
  Field velocityN2;
  Field tractionX;
  Field tractionY;
  Field tractionZ;
  Field stress;
  Field force;
  Field eddyviscosity;
  Field totalviscosity;
  Field utau;
  Field uparallel;
  Field tau;
  Field tauwall; //parallel to flow
  // ESInterface
  Field InterfaceVelocity;
  Field InterfacePressure;
  Field InterfaceStress;
  Field InterfaceDensity;
};

#endif
