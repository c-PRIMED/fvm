// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ICONTAINER_H_
#define _ICONTAINER_H_

#include "misc.h"
#include "RLogInterface.h"
#include "CException.h"

class IContainer
{
public:

  IContainer() {}
  virtual ~IContainer() {}
  
  virtual void zero()  =0;

  // implement operator= in terms of the copyFrom virtual function
  // rather than making it a pure virtual because assignment is always
  // defined in each derived class and we would need to override all versions
  
  IContainer& operator=(const IContainer& oc)
  {
    copyFrom(oc);
    return *this;
  }
  
  virtual shared_ptr<IContainer> newCopy() const = 0;
  virtual shared_ptr<IContainer> newClone() const = 0;

  virtual void copyFrom(const IContainer& oc) = 0;
private:
  IContainer(const IContainer&);
};

#endif
