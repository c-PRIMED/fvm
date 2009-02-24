#ifndef _GEOMFIELDS_H_
#define _GEOMFIELDS_H_

#include "Field.h"

#ifndef SWIG
class StorageSite;
class Matrix;
#endif

struct GeomFields
{
  GeomFields(const string baseName);

  Field coordinate;
  Field area;
  Field areaMag;
  Field volume;
  // this file gets directly included in a swig ineterface definition
  // file hence protect the following
#ifndef SWIG

  
  typedef pair<const StorageSite*, const StorageSite*> SSPair;

  mutable map<SSPair,shared_ptr<Matrix> > _interpolationMatrices;
  
#endif

};

#endif
