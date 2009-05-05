%{
#include "Grid.h"
#include "Array.h"
#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"

  %}

using namespace std;

class Grid
{
 public: 
  Grid();

  ~Grid();

   void Impl(Mesh& mesh, GeomFields& geomFields, FlowFields& flowFields, 
Grid& grid, const string fileBase);	

 void setConnFaceToGrid(Mesh& mesh,  const GeomFields& geomFields, 
			     Grid& grid,  const StorageSite& faces );
	
  const StorageSite& getGrids();

shared_ptr<VecD3Array> computeInterpolatedVelocity(const StorageSite& grids, 
				   Grid& grid,
				   const Mesh& mesh,
				  const GeomFields& geomFields,
	                           const string fileBase, 
				   const StorageSite& faces);

};



 

