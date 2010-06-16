%{
#include "KSearchTree.h"
  %}

class KSearchTree
{
public:

  %extend
  {
    KSearchTree(ArrayBase& pointsBase)
    {
        typedef Vector<double,3> Vec3D;
        typedef Array<Vec3D> Vec3DArray;
        const Vec3DArray& points(dynamic_cast<const Vec3DArray&>(pointsBase));
        return new KSearchTree(points);
    }

    
    void findNeighbors(Vector<double,3>& p, const int k,
                       ArrayBase& neighborsBase)
    {
      typedef Array<int> IntArray;
      IntArray& neighbors(dynamic_cast<IntArray&>(neighborsBase));
      self->findNeighbors(p,k,neighbors);
    }
  }
  
};
