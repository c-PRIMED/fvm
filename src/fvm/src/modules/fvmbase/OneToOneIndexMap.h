#ifndef _ONETOONEINDEXMAP_H_
#define _ONETOONEINDEXMAP_H_

#include "Array.h"

class OneToOneIndexMap 
{
public:
  OneToOneIndexMap(shared_ptr<Array<int> > fromIndices,
                   shared_ptr<Array<int> > toIndices);
  
  
  virtual ~OneToOneIndexMap();

  DEFINE_TYPENAME("OneToOneIndexMap");

  const Array<int>& getFromIndices() const {return *_fromIndices;}
  const Array<int>& getToIndices() const {return *_toIndices;}
  
private:
  shared_ptr<Array<int> > _fromIndices;
  shared_ptr<Array<int> > _toIndices;
};

#endif
