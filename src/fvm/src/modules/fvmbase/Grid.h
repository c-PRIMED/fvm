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
  Grid();

  ~Grid();

  
  
  const StorageSite& getGrids() const {return _grids;}

  const StorageSite& getCells() const {return _gridCells;}
  
  const shared_ptr<Array<VecD3> >& getCoordinates() {return  _coordinates;}
   
  const shared_ptr<Array<VecD3> >& getVelocities() {return  _velocities;}
  
  
  void setCoordinates(const shared_ptr<Array<VecD3> > x) {_coordinates = x;}
 
  void setVelocities(const shared_ptr<Array<VecD3> > x) {_velocities = x;}
 
  const shared_ptr<Array<VecD3> > readCoordinates(const char *file);

  const shared_ptr<Array<VecD3> >  readVelocities(const char *file); 

  void Init (const shared_ptr<Array<VecD3> > coordinates,
	     const shared_ptr<Array<VecD3> > velocities );

  void setandwriteGrids(const string fileBase);

  void setGridMesh(Mesh& mesh);

  void Impl(Mesh& mesh, GeomFields& geomFields, 
	    FlowFields& flowFields, Grid& grid, const string fileBase);

  vector<int> findNeighbors(const VecD3& point, 
			    const Array<VecD3>& grids,
			    const int nGrid);

  const shared_ptr<CRConnectivity>  setConnectivity
                             (const StorageSite& pointSite, 
			      const StorageSite& gridSite,
			      const VecD3Array& points, 
			      const VecD3Array& grids,
			      Grid& grid,
			      const int nGridCells, 
			      const CRConnectivity& cellToGrids );

 void computeInterpolatedVelocity(const StorageSite& grids, 
				   Grid& grid,
				   const Mesh& mesh,
				  const GeomFields& geomFields,
				  const string fileBase);

 vector<int> findNeighborsByCells(const VecD3& point,
		    const Array<VecD3>& grids,
		    const int nCells, 
		    const CRConnectivity& cellToGrids);
 protected:
  StorageSite _grids;
  StorageSite _gridCells;
  shared_ptr<Array<VecD3> > _coordinates;
  shared_ptr<Array<VecD3> > _velocities;


};



#endif
 
