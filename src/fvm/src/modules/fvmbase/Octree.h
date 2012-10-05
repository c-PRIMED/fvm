// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _OCTREE_H_
#define _OCTREE_H_
// -----------------------------------------------------------------------------
// The octree class 
// -----------------------------------------------------------------------------

#include "Array.h"
#include "Vector.h"
#include "Mesh.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include <iostream>
#include <string>
#include <math.h>



class   Octree
{
 public:   
  typedef double T;
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  //this defines a point, which contains cellcentroid coordinate, cellIndex and code
  struct  Point  
  {
        VectorT3      coordinate;    
        int           cellIndex;
        int           code;           // Used during octree generation
  };
// -----------------------------------------------------------------------------
// This defines a cubic bounding volume (center, radius)
// -----------------------------------------------------------------------------

  struct Bounds
  {
        VectorT3        center;         // Center of a cubic bounding volume
        double          radius;         // Radius of a cubic bounding volume
  };




 
  // Construction/Destruction
             Octree();

  virtual   ~Octree();

  // Accessors

  //inline  const   Point * points() const {return _points;}
  //inline  const   unsigned int    pointCount() const {return _pointCount;}

  // Implementation

  //build octree

  shared_ptr<ArrayBase> getArrayPtr(const VectorT3Array&);

  virtual   bool            build(Point *points,
				      const unsigned int count,
                                      const unsigned int threshold,
                                      const unsigned int maximumDepth,
                                      const Bounds &bounds,
                                      const unsigned int currentDepth = 0
				      );

  //calculate bounds
  const   Bounds          calcCubicBounds(const Point * points,
                                                const unsigned int count);

  //report octree
  bool              report(FILE *fp);

  const int               getNode(const double x, const double y, const double z);

  const int               getNode(const VectorT3 coordinate);

  const int               getNode(const VectorT3 coordinate,  double& shortestDistance);

  const double            borderDistance(const VectorT3 coordinate);

  void                    getNodes(const  VectorT3 coordinate,  const double radius, vector<int>& cellList);
  void                    getNodes(const double x, const double y, const double z,
				   const double radius, vector<int>& cellList);
  const int               Naive_getNode(const VectorT3 coordinate, const int count, const Point * points);

  const vector<int>       Naive_getNodes(const VectorT3 coordinate, const int count, const Point * points, const double radius);

  

  void Impl(const Mesh& mesh, const GeomFields& geomFields);

  void Create(const Mesh& mesh, const GeomFields& geomFields, const int faceGroupID);

  //void Impl(const int nCells, shared_ptr<VectorT3Array> cellCentroid); 
//virtual const   bool            traverse(callback proc, void *data) const;

 protected:
  Octree                  *_child[8];         //child node
  unsigned int            _pointCount;        //how many points in this node
  Point                   *_points;           //store the points
  VectorT3                _center;            //node center
  T                       _radius;            //node radius
  unsigned int            _nodeType;          //1=leaf 0=node
  int                     _currentDepth;      //depth or level of node
};
#endif
