// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FIELDLABEL_H_
#define _FIELDLABEL_H_

class FieldLabel
{
public:
  FieldLabel(const string& name) :
    _name(name)
  {}
  
  const string& getName() const {return _name;}
private:
  FieldLabel(const FieldLabel&);
  const string _name;
};

#endif
