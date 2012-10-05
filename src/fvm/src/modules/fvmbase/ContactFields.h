// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CONTACTFIELDS_H_
#define _CONTACTFIELDS_H_
#include "StorageSite.h"
#include "Vector.h"
#include "FloatVarDict.h"
#include "Field.h"
#include "FieldLabel.h"


template<class T>
struct ContactModelConstants : public FloatVarDict<T>
{
  ContactModelConstants()
  {
    this->defineVar("gap", T(3.0e-6));             
    this->defineVar("thickness", T(2.0e-6));    
  }
};

struct ContactFields
{

  ContactFields(const string baseName);

  Field force;
};




#endif 
