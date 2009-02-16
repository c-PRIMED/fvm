#ifndef _OCTREE_H_
#define _OCTREE_H_
// -----------------------------------------------------------------------------
// The octree class 
// -----------------------------------------------------------------------------

#include "Array.h"
#include "Vector.h"

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
  typedef struct   
  {
        VectorT3      coordinate;    
        int           cellIndex;
        int           code;           // Used during octree generation
  } Point;
// -----------------------------------------------------------------------------
// This defines a cubic bounding volume (center, radius)
// -----------------------------------------------------------------------------

  typedef struct
  {
        VectorT3        center;         // Center of a cubic bounding volume
        double          radius;         // Radius of a cubic bounding volume
  } Bounds;




 
  // Construction/Destruction
             Octree();
  virtual   ~Octree();

  // Accessors

  //inline  const   Point * points() const {return _points;}
  //inline  const   unsigned int    pointCount() const {return _pointCount;}

  // Implementation

  //build octree

  virtual const   bool            build(Point *points,
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
  const bool              report(FILE *fp);

  const int               getNode(const double x, const double y, const double z);

  const int               getNode(VectorT3 coordinate);

  const int               getNode(VectorT3 coordinate,  double& shortestDistance);

  const double            borderDistance(VectorT3 coordinate);

  const vector<int>       getNodes(VectorT3 coordinate,  double radius);

  const int               Naive_getNode(VectorT3 coordinate, int count, Point * points);

  const vector<int>       Naive_getNodes(VectorT3 coordinate, int count, Point * points, double radius);

  const bool              MPM_Points_Write(char *file);
  
  const shared_ptr<VectorT3Array>  MPM_Points_Read(char *file);


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
