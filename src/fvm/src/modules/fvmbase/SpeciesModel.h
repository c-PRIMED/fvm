#ifndef _SPECIESMODEL_H_
#define _SPECIESMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "SpeciesFields.h"

#include "Mesh.h"


#include "SpeciesBC.h"

template<class T>
class SpeciesModel : public Model
{
public:

  typedef std::map<int,SpeciesBC<T>*> SpeciesBCMap;
  typedef std::map<int,SpeciesVC<T>*> SpeciesVCMap;
  class Impl;
  
  
  SpeciesModel(const GeomFields& geomFields,
               const MeshList& meshes, 
               const int nSpecies);
  
  virtual ~SpeciesModel();

  virtual void init();
  
  SpeciesFields& getSpeciesFields(const int speciesId);
  SpeciesBCMap& getBCMap(const int speciesId);
  SpeciesVCMap& getVCMap(const int speciesId);
  
  SpeciesBC<T>& getBC(const int id, const int speciesId);

  SpeciesModelOptions<T>& getOptions();

  T getMassFluxIntegral(const Mesh& mesh, const int faceGroupId, const int m);
  bool ResidualLessThanTolerance(const double Tolerance);
  
  void printBCs();

  void updateTime(); 

  void advance(const int niter);
private:
  shared_ptr<Impl> _impl;
};

#endif

