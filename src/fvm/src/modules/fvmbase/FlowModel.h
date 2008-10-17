#ifndef _FLOWMODEL_H_
#define _FLOWMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "FlowFields.h"

#include "Mesh.h"
#include "AMG.h"

#include "FlowBC.h"

template<class T>
class FlowModel : public Model
{
public:

  typedef std::map<int,FlowBC<T>*> FlowBCMap;
  typedef std::map<int,FlowVC<T>*> FlowVCMap;
  
  class Impl;
  
  
  FlowModel(const GeomFields& geomFields,
               FlowFields& thermalFields, const MeshList& meshes);
  
  virtual ~FlowModel();

  virtual void init();
  
  FlowBCMap& getBCMap();
  FlowVCMap& getVCMap();

  FlowModelOptions<T>& getOptions();

  void printBCs();

  void advance(const int niter);

  AMG& getMomentumSolver();
  AMG& getContinuitySolver();

#ifndef USING_ATYPE_TANGENT
  
  void dumpContinuityMatrix(const string fileBase);
  
#endif
  
private:
  shared_ptr<Impl> _impl;
};

#endif

