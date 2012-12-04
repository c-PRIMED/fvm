// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
  typedef Array<int> IntArray;
  typedef Array<VectorT3> VectorT3Array;
      
  MeshMetricsCalculator(GeomFields& geomFields, const MeshList& meshes, bool transient=false);
  virtual ~MeshMetricsCalculator();

  virtual void init();

  void createNodeDisplacement();

  void calculateBoundaryNodeNormal();

  void recalculate();

  void recalculate_deform();

  void computeIBInterpolationMatrices(const StorageSite& particles, const int option=0);

  void computeIBInterpolationMatricesCells();

  void eraseIBInterpolationMatrices(const StorageSite& particles);

  void computeSolidInterpolationMatrices(const StorageSite& particles);

  void computeIBandSolidInterpolationMatrices(const StorageSite& particles);

  void computeGridInterpolationMatrices(const StorageSite& grids, const StorageSite& faces );

  void updateTime();

#ifdef USING_ATYPE_TANGENT
  void setTangentCoords(int meshID, int faceZoneID, int dim);
#endif
private:
  GeomFields& _geomFields;
  Field& _coordField;
  Field& _areaField;
  Field& _areaMagField;
  Field& _volumeField;
  Field& _nodeDisplacement;
  Field& _boundaryNodeNormal;
  bool _transient;

  virtual void calculateNodeCoordinates(const Mesh& mesh);

  virtual void calculateFaceCentroids(const Mesh& mesh);

  virtual void calculateCellCentroids(const Mesh &mesh);
    
  virtual void calculateFaceAreas(const Mesh& mesh);
  
  virtual void calculateFaceAreaMag(const Mesh& mesh);
      
  virtual void calculateCellVolumes(const Mesh& mesh);

  void computeIBInterpolationMatrices(const Mesh& mesh,
                                      const StorageSite& particles,
  				      const int option);
  void computeIBInterpolationMatricesCells(const Mesh& mesh);

  void computeSolidInterpolationMatrices(const Mesh& mesh,
                                         const StorageSite& particles);
  void computeIBandSolidInterpolationMatrices(const Mesh& mesh,
                                      const StorageSite& particles);
  void computeGridInterpolationMatrices(const Mesh& mesh,
  					const StorageSite& grids,
  					const StorageSite& faces  );

};

#endif
