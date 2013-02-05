// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "PhononMacro.h"

PhononMacro::PhononMacro(const std::string baseName) :
  temperature(baseName + ".temperature"),
  deltaT(baseName + ".deltaT"),
  e0(baseName + ".e0"),
  TlResidual(baseName + "TlResidual"),
  TlInjected(baseName + "TlInjected"),
  TlFASCorrection(baseName + "TlFASCorrection"),
  heatFlux(baseName + ".heatFlux"),
  lam(baseName + ".lam"),
  BranchTemperatures(),
  BranchFlux(),
  CellColors(),
  zero(baseName + "zero"),
  one(baseName + "one")
{}

Field& PhononMacro::getModeTemp(int mesh, int mode)
{
  FieldVectorMap::iterator it;
  it=BranchTemperatures.find(mesh);
  if(it!=BranchTemperatures.end())
    {
      FieldVector& FieldVec=*(it->second);
      if(mode<FieldVec.size())
	return *FieldVec[mode];
      else
	throw CException("Don't have Field for that mode!");
    }
  throw CException("Don't have FieldVector for that mesh!");
}

Field& PhononMacro::getModeFlux(int mesh, int mode)
{
  FieldVectorMap::iterator it;
  it=BranchFlux.find(mesh);
  if(it!=BranchFlux.end())
    {
      FieldVector& FieldVec=*(it->second);
      if(mode<FieldVec.size())
	return *FieldVec[mode];
      else
	throw CException("Don't have Field for that mode!");
    }
  throw CException("Don't have FieldVector for that mesh!");
}

Field& PhononMacro::getColorField(int level)
{return *CellColors[level];}

Field& PhononMacro::getPlotColorField(int level)
{return *plottingCellColors[level];}
