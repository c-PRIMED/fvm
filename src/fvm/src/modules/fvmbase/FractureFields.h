// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FRACTUREFIELDS_H_
#define _FRACTUREFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct FractureFields
{
  FractureFields(const string baseName);

  Field phasefieldvalue;
  Field phasefieldvalueN1;
  Field phasefieldvalueN2;
  //Field specificHeat;
  Field phasefieldFlux;
  Field phasefieldGradient;
  Field conductivity;
  Field source;
  Field sourcecoef;
  //Field convectionFlux;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
};

#endif
