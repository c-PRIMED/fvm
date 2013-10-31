// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _VACANCYMODEL_H_
#define _VACANCYMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "VacancyFields.h"

#include "Mesh.h"


#include "VacancyBC.h"

template<class T>
class VacancyModel : public Model
{
public:

  typedef std::map<int,VacancyBC<T>*> VacancyBCMap;
  typedef std::map<int,VacancyVC<T>*> VacancyVCMap;
  class Impl;
  
  
  VacancyModel(const GeomFields& geomFields,
               VacancyFields& vacancyFields, const MeshList& meshes);
  
  virtual ~VacancyModel();

  virtual void init();
  
  VacancyBCMap& getBCMap();
  VacancyVCMap& getVCMap();
  
  VacancyBC<T>& getBC(const int id);

  VacancyModelOptions<T>& getOptions();
  
  void printBCs();

  void advance(const int niter);

  void updateTime();
private:
  shared_ptr<Impl> _impl;
};

#endif

