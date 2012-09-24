#ifndef _MACROPARAMETERS_H_
#define _MACROPARAMETERS_H_


#include "Field.h"
#include "FlowFields.h"

struct MacroFields //public Flowfields
{
  MacroFields(const string baseName);//: Flowfields(baseName){
  Field velocity;
  Field velocityResidual;
  Field velocityInjected;
  Field velocityFASCorrection;
  Field pressure;
  Field viscosity;
  Field density;
  Field temperature;
  Field collisionFrequency;
  Field Txx;
  Field Tyy;
  Field Tzz;
  Field Txy;
  Field Tyz;
  Field Tzx;
  Field coeff;
  Field coeffg;
  Field Entropy;
  Field EntropyGenRate;
  Field EntropyGenRate_Collisional;
  Field force;
  Field Stress;
  Field heatFlux;
  //Field distribution;
  //Field distributionGradient;
  
  Field InterfaceVelocity;
  Field InterfacePressure;
  Field InterfaceStress;
  Field InterfaceDensity;

  //Field M300; //M300,M120,M102 cx^3,    cx*cy^2,cx*cz^2
  //Field M030; //M210 M030 M012 cy*cx^2,   cy^3, cy*cz^2
  //Field M003; //M201,M021,M003 cz*cx^2, cz*cy^2, cz^3
  Field Knq;  //0 for x-dir,1 for y-dir, 2 for z-dir variation
  
};

#endif
