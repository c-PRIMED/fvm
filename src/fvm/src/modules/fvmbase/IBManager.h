#ifndef _IBMANAGER_H_
#define _IBMANAGER_H_

#include "Mesh.h"
#include "GeomFields.h"
#include <set>

class KSearchTree;
class AABB;

/**
 * structure that keeps information about nearest cell for a solid
 * boundary face.
 * 
 */

struct
NearestCell
{
  NearestCell():
    mesh(0),
    cell(-1),
    distanceSquared(0)
  {}

  const Mesh* mesh;
  int cell;
  double distanceSquared;
  set<int> neighbors;
};

struct   
doubleIntStruct{
   double VALUE;
   int TAG;
};



class IBManager
{
public:
  typedef Vector<double,3> Vec3D;
  typedef Array<Vec3D> Vec3DArray;
  typedef Array<int> IntArray;

  IBManager(GeomFields& geomFields,
            Mesh& solidBoundaryMesh,
            const MeshList& fluidMeshes);

  void update();

  int fluidNeighborsPerIBFace;
  int fluidNeighborsPerSolidFace;
  int solidNeighborsPerIBFace;
  
private:

  void markIntersections(Mesh& fluidMesh, AABB& sMeshesAABB);
  int markFluid(Mesh& fluidMesh);
  int markSolid(Mesh& fluidMesh);
  void markIBTypePlus(Mesh& fluidMesh);
  void createIBFaces(Mesh& fluidMesh);
  void createIBInterpolationStencil(Mesh& mesh,
                                    KSearchTree& fluidCellsTree,
                                    KSearchTree& solidFacesTree);
  void
  findNearestCellForSolidFaces(Mesh& mesh,
                               KSearchTree& fluidCellsTree,
                               vector<NearestCell>& nearest);
  
  void createSolidInterpolationStencil(Mesh& mesh,
                                       vector<NearestCell>& nearest);                                       
  GeomFields& _geomFields;
  Mesh& _solidBoundaryMesh;
  const MeshList _fluidMeshes;
};
#endif
