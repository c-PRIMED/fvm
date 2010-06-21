#ifndef _FLOWMODEL_H_
#define _FLOWMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "FlowFields.h"

#include "Mesh.h"
#include "LinearSolver.h"

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
  
  virtual map<string,shared_ptr<ArrayBase> >&
  getPersistenceData();

  virtual void restart();
  
  FlowBCMap& getBCMap();
  FlowVCMap& getVCMap();

  FlowModelOptions<T>& getOptions();

  void printBCs();

  // do the specified number of iterations, return true if converged 
  bool advance(const int niter);
  bool advanceCoupled(const int niter);

  void updateTime();

  Vector<T,3> getPressureIntegral(const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getPVIntegral(const Field& velCoeff, const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getMomentumFluxIntegral(const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getMomentumDerivativeIntegral(const Mesh& mesh);


  Vector<T,3> getPressureIntegralonIBFaces(const Mesh& mesh);
  Vector<T,3> getMomentumFluxIntegralonIBFaces(const Mesh& mesh);


  boost::shared_ptr<ArrayBase> getStressTensor(const Mesh& mesh, const ArrayBase& cellIds);

  
  void printPressureIntegrals();
  void printMomentumFluxIntegrals();
  void printMassFluxIntegrals();
  
  void computeIBFaceVelocity(const StorageSite& particles);
  void computeSolidSurfaceStress(const StorageSite& particles);
  void computeIBandSolidVelocity(const StorageSite& particles);
  //LinearSolver& getMomentumSolver();
  //LinearSolver& getContinuitySolver();

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
  
  void dumpContinuityMatrix(const string fileBase);
  
#endif
  
private:
  shared_ptr<Impl> _impl;
};

#endif

