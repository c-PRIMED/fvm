// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GLOBALFIELDS_H_
#define _GLOBALFIELDS_H_

#include "Field.h"

#include "FieldLabel.h"

class GlobalFields
{
public:
  typedef map<const FieldLabel*, shared_ptr<Field> > FieldMap;
  static Field& createField(const FieldLabel& fieldLabel);
  static const Field& getField(const FieldLabel& fieldLabel);

private:

  static FieldMap& getFieldMap();
  static FieldMap* _fieldMap;
};

#define getGlobalField GlobalFields::getField
#define createGlobalField GlobalFields::createField

#endif
