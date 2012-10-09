// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECTRICFIELDS_H_
#define _ELECTRICFIELDS_H_


#include "Field.h"
#include "FieldLabel.h"

struct ElectricFields
{
  ElectricFields(const string baseName);

  //Fields in electrastatics

  Field potential;
  Field potential_flux;              /* this is only stored on boundary faces for the purpose */
                                     /* of post processing and applying generic bc*/
  Field potential_gradient;
  Field electric_field;              //potentialGradient vector; 
  Field dielectric_constant;         //permittivity;
  Field total_charge;                 //source term in Poisson equation
                                     //which is the sum of charge[0] and charge[1]                           
  

  //Fields in charge transport

  Field conduction_band;
  Field valence_band;
  Field electron_totaltraps;          //number of electron traps at each cell
  Field free_electron_capture_cross; //free electron capture cross section
  Field transmission;
  Field electron_velocity;
  Field charge;
  Field chargeFlux;
  Field diffusivity;
  Field convectionFlux;
  Field chargeGradient;
  Field chargeN1;
  Field chargeN2;
  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density
 
  Field force;
  //Field oneD_column;              //the one dimensional columns used in dielectric chargine 1D model
  
  //Fields for coupling to species model
  Field speciesConcentration;
  Field lnSpeciesConcentration;
  Field lnSpeciesConcentrationGradient;

};

#endif
