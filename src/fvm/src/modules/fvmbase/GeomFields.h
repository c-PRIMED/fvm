#ifndef _GEOMFIELDS_H_
#define _GEOMFIELDS_H_

#include "Field.h"

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

  
#endif

};

#endif
