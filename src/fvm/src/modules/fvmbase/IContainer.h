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
  virtual IContainer& operator=(const IContainer& oc)  =0;
  
  virtual shared_ptr<IContainer> newCopy() const = 0;
  virtual shared_ptr<IContainer> newClone() const = 0;

};

#endif
