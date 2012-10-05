// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GRID_H_
#define _GRID_H_

#include "Array.h"
#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "Mesh.h"
#include "GeomFields.h"
#include "FlowFields.h"
#include "MPM_Particles.h"
#include "CRMatrixTranspose.h"



typedef Vector<double,3> VecD3;
  
typedef Array<VecD3> VecD3Array;

class Grid
{
public: 
  Grid(GeomFields& geomFields, FlowFields& flowFields, string coordFileName,
       string velocityFileName);

  ~Grid();

  const shared_ptr<Array<VecD3> >& getCoordinates() {return  _coordinates;}
   
  const shared_ptr<Array<VecD3> >& getVelocities() {return  _velocities;}
  
  const StorageSite& getNodes() {return _nodes;}
  const StorageSite& getCells() {return _cells;}
  

  void setandwriteGrids(const string fileBase);

  void createCellToNodeConnectivity();

  void setConnFaceToGrid(Mesh& mesh, const StorageSite& faces);

  vector<int> findNeighbors(const VecD3& point);

  const shared_ptr<CRConnectivity>
  createConnectivity(const StorageSite& pointSite, 
                     const VecD3Array& points);

  shared_ptr<ArrayBase>
  computeInterpolatedVelocity(const StorageSite& faces);
  
  vector<int>
  findNeighborsByCells(const VecD3& point);
  
protected:
  GeomFields& _geomFields;
  FlowFields& _flowFields;
  StorageSite _nodes;
  StorageSite _cells;
  shared_ptr<Array<VecD3> > _coordinates;
  shared_ptr<Array<VecD3> > _velocities;
  shared_ptr<CRConnectivity> _cellNodes;

};



#endif
 
