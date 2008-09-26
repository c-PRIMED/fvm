#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "ThermalFields.h"

#include "Mesh.h"


#include "ThermalBC.h"

template<class T>
class ThermalModel : public Model
{
public:

  typedef std::map<int,ThermalBC<T>*> ThermalBCMap;
  class Impl;
  
  
  ThermalModel(const GeomFields& geomFields,
               ThermalFields& thermalFields, const MeshList& meshes);
  
  virtual ~ThermalModel();

  virtual void init();
  
  ThermalBCMap& getBCMap();
  ThermalBC<T>& getBC(const int id);

  ThermalModelOptions<T>& getOptions();

  void printBCs();

  void advance(const int niter);
private:
  shared_ptr<Impl> _impl;
};

#endif

