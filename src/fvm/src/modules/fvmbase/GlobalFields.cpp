// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "GlobalFields.h"

GlobalFields::FieldMap *GlobalFields::_fieldMap = 0;

Field&
GlobalFields::createField(const FieldLabel& fieldLabel)
{
  FieldMap& fMap = getFieldMap();
  FieldMap::iterator p = fMap.find(&fieldLabel);
  if (p == fMap.end())
  {
      shared_ptr<Field> f(new Field(fieldLabel.getName()));
      fMap[&fieldLabel] = f;
      return *f;
  }
  else
    return *p->second; 
}

const Field&
GlobalFields::getField(const FieldLabel& fieldLabel)
{
  FieldMap& fMap = getFieldMap();
  FieldMap::iterator p = fMap.find(&fieldLabel);
  if (p != fMap.end()) return *p->second;
  throw CException("No such field " + fieldLabel.getName());
}

GlobalFields::FieldMap&
GlobalFields::getFieldMap()
{
  if (!_fieldMap)
  {
      _fieldMap = new FieldMap();
  }
  return *_fieldMap;
}
