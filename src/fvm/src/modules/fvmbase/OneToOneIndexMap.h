#ifndef _ONETOONEINDEXMAP_H_
#define _ONETOONEINDEXMAP_H_

#include "Array.h"

class OneToOneIndexMap 
{
public:
  OneToOneIndexMap(const Array<int>* fromIndices,
                   const Array<int>* toIndices);
  
  
  virtual ~OneToOneIndexMap();

  DEFINE_TYPENAME("OneToOneIndexMap");

  const Array<int>& fromIndices;
  const Array<int>& toIndices;
};

#endif
