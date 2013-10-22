// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BATTERYMODEL_H_
#define _BATTERYMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "BatteryFields.h"

#include "Mesh.h"


#include "BatteryBC.h"

template<class T>
class BatteryModel : public Model
{
public:

  typedef std::map<int,BatterySpeciesBC<T>*> BatterySpeciesBCMap;
  typedef std::map<int,BatterySpeciesVC<T>*> BatterySpeciesVCMap;
  typedef std::map<int,BatteryPotentialBC<T>*> BatteryPotentialBCMap;
  typedef std::map<int,BatteryPotentialVC<T>*> BatteryPotentialVCMap;
  typedef std::map<int,BatteryThermalBC<T>*> BatteryThermalBCMap;
  typedef std::map<int,BatteryThermalVC<T>*> BatteryThermalVCMap;

  class Impl;
  
  
  BatteryModel(GeomFields& geomFields,
	       const MeshList& realMeshes, 
               const MeshList& meshes, 
               const int nSpecies);
  
  virtual ~BatteryModel();

  virtual void init();
  
  BatterySpeciesFields& getBatterySpeciesFields(const int speciesId);
  BatteryModelFields& getBatteryModelFields();
  BatterySpeciesBCMap& getSpeciesBCMap(const int speciesId);
  BatteryPotentialBCMap& getPotentialBCMap();
  BatteryThermalBCMap& getThermalBCMap();
  BatterySpeciesVCMap& getSpeciesVCMap(const int speciesId);
  BatteryPotentialVCMap& getPotentialVCMap();
  BatteryThermalVCMap& getThermalVCMap();
  
  //SpeciesBC<T>& getBC(const int id, const int speciesId);

  BatteryModelOptions<T>& getOptions();

  T getMassFluxIntegral(const Mesh& mesh, const int faceGroupId, const int m);
  T getPotentialFluxIntegral(const Mesh& mesh, const int faceGroupId);
  T getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId);
  T getAverageConcentration(const Mesh& mesh, const int m);
  T getFaceGroupArea(const Mesh& mesh, const int fgID);
  T getFaceGroupVoltage(const Mesh& mesh, const int fgID);
  T getMeshVolume(const Mesh& mesh);
  T getSpeciesResidual(const int speciesId);
  T getPotentialResidual();
  T getThermalResidual();
  T getPCResidual(const int v);
  
  //void printBCs();

  void updateTime(); 
  void recoverLastTimestep(); 
  void advanceSpecies(const int niter);
  void advancePotential(const int niter);
  void advanceThermal(const int niter);
  void advanceCoupled(const int niter);

  void copySeparateToCoupled();
  void copyCoupledToSeparate();
private:
  shared_ptr<Impl> _impl;
};

#endif

