%{
#include "IBManager.h"
  %}

class IBManager
{
public:

  IBManager(GeomFields& geomFields,
            Mesh& solidBoundaryMesh,
            const MeshList& fluidMeshes);

  void update();

  int fluidNeighborsPerIBFace;
  int solidNeighborsPerIBFace;
  int fluidNeighborsPerSolidFace;
  int IBNeighborsPerSolidFace;
};
