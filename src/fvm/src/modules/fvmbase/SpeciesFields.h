#ifndef _SPECIESFIELDS_H_
#define _SPECIESFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct SpeciesFields
{
  SpeciesFields(const string baseName);

  Field massFraction;
  Field massFlux;
  Field diffusivity;
  Field source;
  Field convectionFlux;

  Field massFractionN1;
  Field massFractionN2;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
};

struct SpeciesModelFields
{
  SpeciesModelFields(const string baseName);

  Field speciesGradient;
};

#endif
