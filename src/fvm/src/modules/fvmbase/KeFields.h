// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KEFIELDS_H_
#define _KEFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct KeFields
{
  KeFields(const string baseName);

  Field energy;
  Field energyGradient;
  Field dissipation;
  Field dissipationGradient;
  Field c1;  //eddy viscosi/sigmak
  Field c2;  //eddy viscosi/sigmae
  Field sigmak;
  Field sigmae;
  Field cmu;
  //Field convectionFlux;
  Field kFlux;
  Field eFlux;
  Field energyN1;
  Field energyN2;
  //Field density;
  Field dissipationN1;
  Field dissipationN2;
  Field sourcek;   // S of energy
  Field sourced;  // S of dissipation
  Field sourcec;  //Sc after linearization
  Field sourcep;  //Sp after linearization
};

#endif




