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
  BandTemperatures(),
  BandRelEnergy(),
  BandFlux(),
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

Field& PhononMacro::getBandTemp(int mesh, int band)
{
  FieldVectorMap::iterator it;
  it=BandTemperatures.find(mesh);
  if(it!=BranchTemperatures.end())
    {
      FieldVector& FieldVec=*(it->second);
      if(band<FieldVec.size())
	return *FieldVec[band];
      else
	throw CException("Band Temp: Don't have Field for that mode!");
    }
  throw CException("Band Temp: Don't have FieldVector for that mesh!");
}

Field& PhononMacro::getBandRelEnergy(int mesh, int band)
{
  FieldVectorMap::iterator it;
  it=BandRelEnergy.find(mesh);
  if(it!=BandRelEnergy.end())
    {
      FieldVector& FieldVec=*(it->second);
      if(band<FieldVec.size())
	return *FieldVec[band];
      else
	throw CException("Band Relative Energy: Don't have Field for that mode!");
    }
  throw CException("Band Relative Energy: Don't have FieldVector for that mesh!");
}

Field& PhononMacro::getBandFlux(int mesh, int band)
{
  FieldVectorMap::iterator it;
  it=BandFlux.find(mesh);
  if(it!=BandFlux.end())
    {
      FieldVector& FieldVec=*(it->second);
      if(band<FieldVec.size())
	return *FieldVec[band];
      else
	throw CException("Band Flux: Don't have Field for that mode!");
    }
  throw CException("Band Flux: Don't have FieldVector for that mesh!");
}

Field& PhononMacro::getColorField(int level)
{return *CellColors[level];}

Field& PhononMacro::getPlotColorField(int level)
{return *plottingCellColors[level];}
