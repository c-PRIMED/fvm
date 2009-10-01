#ifndef _ELECTRONICFIELDS_H_
#define _ELECTRONICFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct ElectronicFields
{
  ElectronicFields(const string baseName);
  //potential field
  Field potential;
  Field potentialFlux;
  Field potentialGradient;           //electronic field vector
  Field permittivity;
  Field source;
  //charge field
  Field charge;
  Field chargeFlux;
  Field diffusivity;
  Field convectionFlux;
  Field chargeGradient;
  Field chargeN1;
  Field chargeN2;
  //tunneling charge
  Field tunnelingCharge;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
 
};

#endif
