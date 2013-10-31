// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _VACANCYFIELDS_H_
#define _VACANCYFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct VacancyFields
{
  VacancyFields(const string baseName);

  Field concentration;
  Field concentrationN1;
  Field concentrationN2;
  Field vacaFlux;
  Field concentrationGradient;
  Field diffusioncoefficient;
 
  Field one;                      //used to fill in density
};

#endif
