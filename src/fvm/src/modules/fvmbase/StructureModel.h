#ifndef _STRUCTUREMODEL_H_
#define _STRUCTUREMODEL_H_

#include "Model.h"

#include "GeomFields.h"
#include "StructureFields.h"

#include "Mesh.h"
#include "LinearSolver.h"

#include "StructureBC.h"

template<class T>
class StructureModel : public Model
{
public:

  typedef std::map<int,StructureBC<T>*> StructureBCMap;
  typedef std::map<int,StructureVC<T>*> StructureVCMap;
  
  class Impl;
  
  
  StructureModel(const GeomFields& geomFields,
               StructureFields& structureFields, const MeshList& meshes);
  
  virtual ~StructureModel();

  virtual void init();
  
  StructureBCMap& getBCMap();
  StructureVCMap& getVCMap();

  StructureModelOptions<T>& getOptions();

  void printBCs();

  // do the specified number of iterations, return true if converged 
  bool advance(const int niter);
  //  bool advanceCoupled(const int niter);

  void updateTime();
  void updateForceOnBoundary(const StorageSite& faceSite, const ArrayBase& bforceA, const map<int,int>& commonFacesMap, 
                             ArrayBase& fxA, ArrayBase& fyA, ArrayBase& fzA);

  //  Vector<T,3> getPressureIntegral(const Mesh& mesh, const int faceGroupID);
  //  Vector<T,3> getPVIntegral(const Field& velCoeff, const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getDeformationFluxIntegral(const Mesh& mesh, const int faceGroupID);
  Vector<T,3> getDeformationDerivativeIntegral(const Mesh& mesh);
  void getTraction(const Mesh& mesh);

  //  Vector<T,3> getPressureIntegralonIBFaces(const Mesh& mesh);
  //  Vector<T,3> getMomentumFluxIntegralonIBFaces(const Mesh& mesh);


  boost::shared_ptr<ArrayBase> getStressTensor(const Mesh& mesh, const ArrayBase& cellIds);

  
  //  void printPressureIntegrals();
  void printDeformationFluxIntegrals();
  //  void printMassFluxIntegrals();
  
  //  void computeIBFaceVelocity(const StorageSite& particles);
  //  void computeIBandSolidVelocity(const StorageSite& particles);
  //LinearSolver& getMomentumSolver();
  //LinearSolver& getContinuitySolver();

private:
  shared_ptr<Impl> _impl;
};

#endif

