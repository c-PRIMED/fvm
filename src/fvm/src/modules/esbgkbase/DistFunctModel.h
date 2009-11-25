#ifndef _DISTFUNCTMODEL_H_
#define _DISTFUNCTMODEL_H_

#include "Model.h"
#include "Array.h"
#include "Vector.h"
#include "quadrature.h"

template<class T>
class DistFunctModel : public Model
{
public:
  typedef Vector<T,3> VectorT3; 
  
  
  DistFunctModel(const GeomFields& geomFields,
               FlowFields& thermalFields, const MeshList& meshes);
  
  virtual ~FlowModel();

  virtual void init();
  
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
  void computeIBandSolidVelocity(const StorageSite& particles);
  //LinearSolver& getMomentumSolver();
  //LinearSolver& getContinuitySolver();

 
private:
  shared_ptr<Impl> _impl;
};

#endif

