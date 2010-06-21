#ifndef _IBMANAGER_H_
#define _IBMANAGER_H_

#include "Mesh.h"
#include "GeomFields.h"

class KSearchTree;
class AABB;

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

  void markIB(Mesh& fluidMesh, AABB& sMeshesAABB);
  void createIBInterpolationStencil(Mesh& mesh,
                                    KSearchTree& fluidCellsTree,
                                    KSearchTree& solidFacesTree);
  void createSolidInterpolationStencil(Mesh& mesh,
                                       KSearchTree& fluidCellsTree);
                                       
  GeomFields& _geomFields;
  Mesh& _solidBoundaryMesh;
  const MeshList _fluidMeshes;
};
#endif
