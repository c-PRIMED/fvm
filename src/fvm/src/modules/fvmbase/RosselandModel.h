// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ROSSELANDMODEL_H_
#define _ROSSELANDMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "FlowFields.h"
#include "ThermalFields.h"
#include "Mesh.h"
#include "LinearSolver.h"


#include "misc.h"
#include "FloatVarDict.h"


//thermalFields =  ThermalFields('therm')

template<class T>
struct RosselandVC : public FloatVarDict<T>
{
  RosselandVC()
  {
    this->defineVar("temperature",T(300.00));
    //this->defineVar("temperature",_thermalFields.temperature);

  }
};


template<class T>
class RosselandModel : public Model
{
public:

  typedef std::map<int,RosselandVC<T>*> VCMap;

  class Impl;


  RosselandModel(const GeomFields& geomFields,
               ThermalFields& thermalFields, const MeshList& meshes);

  virtual ~RosselandModel();

  virtual void init();

  VCMap& getVCMap();

  // do the specified number of iterations, return true if converged
  bool advance(const int niter);

private:
  shared_ptr<Impl> _impl;
};

#endif


