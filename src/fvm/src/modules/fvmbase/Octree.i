%{

#include "Array.h"
//#include "Vector.h"

#include <iostream>
#include <string>
#include <math.h>
#include "ArrayBase.h"
#include "Octree.h"
  %}


using namespace std;

class   Octree
{
public:
  //typedef double T;
  //typedef Array<T> TArray;
  //typedef Vector<T,3> VectorT3;
  //typedef Array<VectorT3> VectorT3Array;	
 

             Octree();
  virtual   ~Octree();


  struct Point   
  {
    //      VecD3      coordinate;    
        int           cellIndex;
        int           code;           // Used during octree generation
  };

  //typedef point_struct Point;

  struct Bound
  {
    //  VecD3        center;         // Center of a cubic bounding volume
        double          radius;         // Radius of a cubic bounding volume
  };

   void Impl(const Mesh& mesh, const GeomFields& geomFields);
 //  void Impl(const int nCells, shared_ptr<VecD3Array> cellCentroid);

   void Create(const Mesh& mesh, const GeomFields& geomFields, const int faceGroupID);

   const int               getNode(const double x, const double y, const double z);

   const int               getNode(const VecD3 coordinate);

   const int               getNode(const VecD3 coordinate,  double& shortestDistance);

   void getNodes(const VecD3 coordinate,  const double radius, vector<int>& cellList);
};


 


