#include "OneToOneIndexMap.h"

OneToOneIndexMap::OneToOneIndexMap(const Array<int>* fromIndices,
                                   const Array<int>* toIndices) :
  fromIndices(*fromIndices),
  toIndices(*toIndices)
{
  logCtor();
}

  
OneToOneIndexMap::~OneToOneIndexMap()
{
  logDtor();
}

