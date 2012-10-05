// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

#include "PyCreatable.h"

class StorageSite;

#include "base_public.h"

class  BASE_PUBLIC Connectivity : public PyCreatable
{
public:
  Connectivity(const Args& args);
  Connectivity(State state,
               const StorageSite& rowSite, const StorageSite& colSite);
  
  virtual ~Connectivity();

  DECLARE_HT("Connectivity");
  DECLARE_METHOD(createOffset);
  
protected:
  mutable StorageSite const *_rowSite;
  mutable StorageSite const *_colSite;
  mutable int _rowDim;
  mutable int _colDim;
};

#endif
