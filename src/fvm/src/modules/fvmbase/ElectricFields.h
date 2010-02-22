#ifndef _ElectricFIELDS_H_
#define _ElectricFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct ElectricFields
{
  ElectricFields(const string baseName);
  //electric field
  Field potential;
  Field potential_flux;              /*this is only stored on boundary faces for the purpose of 
				     post processing and applying generic bc*/
  Field electric_field;              //potentialGradient vector; 
  Field dielectric_constant;         //permittivity;
  Field total_charge;                //source term in Poisson equation
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
