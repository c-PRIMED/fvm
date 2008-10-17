#ifndef _UMESHMETRICSCALCULATOR_H_
#define _UMESHMETRICSCALCULATOR_H_

#include "atype.h"
#include "Model.h"

#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "GeomFields.h"

#include "Mesh.h"

template<class T>
class MeshMetricsCalculator : public Model
{
public:

  typedef Vector<T,3> VectorT3;
  typedef Array<T> TArray;
  typedef Array<VectorT3> VectorT3Array;
      
  MeshMetricsCalculator(GeomFields& geomFields, const MeshList& meshes);
  virtual ~MeshMetricsCalculator();

  virtual void init();

  DEFINE_TYPENAME("MeshMetricsCalculator<"+NumTypeTraits<T>::getTypeName()+">");

#ifdef USING_ATYPE_TANGENT
  void setTangentCoords(int meshID, int faceZoneID, int dim);
#endif
private:
  Field& _coordField;
  Field& _areaField;
  Field& _areaMagField;
  Field& _volumeField;

  virtual void calculateNodeCoordinates(const Mesh& mesh);

  virtual void calculateFaceCentroids(const Mesh& mesh);

  virtual void calculateCellCentroids(const Mesh &mesh);
    
  virtual void calculateFaceAreas(const Mesh& mesh);
  
  virtual void calculateFaceAreaMag(const Mesh& mesh);
      
  virtual void calculateCellVolumes(const Mesh& mesh);

};

#endif
