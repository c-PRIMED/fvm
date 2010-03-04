#ifndef _STRUCTUREFIELDS_H_
#define _STRUCTUREFIELDS_H_


#include "Field.h"

struct StructureFields
{
  StructureFields(const string baseName);

  Field deformation;
  Field deformationGradient;
  Field deformationFlux;
  Field eta;
  Field eta1;
  Field density;
  Field deformationN1;
  Field deformationN2;
};

#endif
