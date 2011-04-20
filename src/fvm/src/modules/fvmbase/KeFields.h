#ifndef _KEFIELDS_H_
#define _KEFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct KeFields
{
  KeFields(const string baseName);

  Field energy;
  Field energyGradient;
  Field dissipation;
  Field dissipationGradient;
  Field c1;  //eddy viscosi/sigmak
  Field c2;  //eddy viscosi/sigmae
  Field sigmak;
  Field sigmae;
  Field cmu;
  //Field convectionFlux;
  Field kFlux;
  Field eFlux;
  Field energyN1;
  Field energyN2;
  //Field density;
  Field dissipationN1;
  Field dissipationN2;
};

#endif




