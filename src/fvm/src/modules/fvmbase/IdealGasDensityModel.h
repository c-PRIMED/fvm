// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _IDEALGASDENSITYMODEL_H_
#define _IDEALGASDENSITYMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "FlowFields.h"

#include "Mesh.h"
#include "LinearSolver.h"


#include "misc.h"
#include "FloatVarDict.h"



template<class T>
struct IdealGasVC : public FloatVarDict<T>
{
  IdealGasVC()
  {
    this->defineVar("pressure",T(0.0));
    this->defineVar("operatingPressure",T(101325.0));
    this->defineVar("temperature",T(300.0));
    this->defineVar("molecularWeight",T(28.966));
    this->defineVar("urf",T(1.0));
    
  }
};


template<class T>
class IdealGasDensityModel : public Model
{
public:

  typedef std::map<int,IdealGasVC<T>*> VCMap;
  
  class Impl;
  
  
  IdealGasDensityModel(const GeomFields& geomFields,
               FlowFields& thermalFields, const MeshList& meshes);
  
  virtual ~IdealGasDensityModel();

  virtual void init();
  
  VCMap& getVCMap();

  // do the specified number of iterations, return true if converged 
  bool advance(const int niter);

private:
  shared_ptr<Impl> _impl;
};


#endif
