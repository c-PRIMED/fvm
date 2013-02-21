// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PHONONMACRO_H_
#define _PHONONMACRO_H_

#include "Field.h"

using namespace std;

struct PhononMacro
{

  typedef shared_ptr<Field> FieldPtr;
  typedef vector<FieldPtr> FieldVector;
  typedef map<int,FieldVector*> FieldVectorMap;

  PhononMacro(const string baseName);

  Field temperature;
  Field deltaT;
  Field e0;
  Field TlResidual;
  Field TlInjected;
  Field TlFASCorrection;
  Field heatFlux;
  Field lam;
  FieldVectorMap BranchTemperatures;
  FieldVectorMap BranchFlux;
  FieldVectorMap BandTemperatures;
  FieldVectorMap BandRelEnergy;
  FieldVectorMap BandFlux;
  FieldVector CellColors;
  FieldVector plottingCellColors;

  Field zero;                     //used to fill in continuityResidual
  Field one;                      //used to fill in density

  Field& getModeTemp(int mesh, int mode);
  Field& getModeFlux(int mesh, int mode);
  Field& getBandTemp(int mesh, int band);
  Field& getBandRelEnergy(int mesh, int band);
  Field& getBandFlux(int mesh, int band);
  Field& getColorField(int level);
  Field& getPlotColorField(int level);
};

#endif
