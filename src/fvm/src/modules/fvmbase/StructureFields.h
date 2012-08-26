#ifndef _STRUCTUREFIELDS_H_
#define _STRUCTUREFIELDS_H_


#include "Field.h"

struct StructureFields
{
  StructureFields(const string baseName);

  Field deformation;
  Field elasticDeformation;
  Field deformationGradient;
  Field deformationFlux;
  Field eta;
  Field eta1;
  Field alpha;
  Field density;
  Field deformationN1;
  Field deformationN2;
  Field deformationN3;
  Field tractionX;
  Field tractionY;
  Field tractionZ;
  Field strainX;
  Field strainY;
  Field strainZ;
  Field plasticDiagStrain;
  Field devStress;
  Field VMStress;
  Field plasticStrain;
  Field temperature;
  Field bodyForce;
  Field volume0;
};

#endif
