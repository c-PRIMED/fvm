// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "OneToOneIndexMap.h"

OneToOneIndexMap::OneToOneIndexMap(shared_ptr<Array<int> > fromIndices,
                                   shared_ptr<Array<int> > toIndices) :
  _fromIndices(fromIndices),
  _toIndices(toIndices)
{
  logCtor();
}

  
OneToOneIndexMap::~OneToOneIndexMap()
{
  logDtor();
}

