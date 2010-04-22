#ifndef _DIELECTRICCOLUMN_H_
#define _DIELECTRICCOLUMN_H_

#include "Vector.h"
#include "Array.h"
#include "ElectricUtilityFunctions.h"
#include "StorageSite.h"
#include "CRMatrixTranspose.h"
#include "PhysicsConstant.h"
#include "CException.h"


template <class T>
class DielectricOneDimColumn
{
 public:
  typedef Vector<T, 3> VectorT3;
  typedef Vector<double, 3> VectorD3;
  
  DielectricOneDimColumn(const StorageSite& cells, 
			 const int nGrid):
    _startPoint(),
    _endPoint(),
    _normalizedDirection(),
    _radius(),
    _grids(nGrid),
    _cells(cells),
    _gridList(nGrid),
    _cellList(0),
    _conduction_band(nGrid),
    _transmission(nGrid),
    _cellGrids(),
    _gridCells(), 
    _cellToGrid(), 
    _gridToCell()
      {   }

  ~DielectricOneDimColumn() {}

   
  void setStartPoint(VectorD3 point) {_startPoint = point;}
  void setEndPoint(VectorD3 point) {_endPoint = point;}
  void setNormalizedDirection(VectorD3 direction) {_normalizedDirection = direction;}
  void setCellList (const int i, const int j) {_cellList[i] = j;}
  void setRadius (const T rr)     { _radius = rr;}
  //void setCells(const StorageSite& cells) { _cells = cells;}

  const VectorD3 getStartPoint() const { return _startPoint;}
  const VectorD3 getEndPoint() const {return _endPoint; }
  const VectorD3 getNormalizedDirection() const  {return _normalizedDirection;}
  const Array<int>& getCellList() const {return _cellList;}

  VectorD3 getStartPoint() { return _startPoint;}
  VectorD3 getEndPoint() {return _endPoint; }
  VectorD3 getNormalizedDirection() {return _normalizedDirection;}
  Array<int>& getCellList() {return _cellList;}

  void discretizeCenterLine(){
    const int nGrid = _grids.getCount();
    double length = mag(_startPoint-_endPoint);
    double delta = length/(nGrid-1);
    for(int n=0; n<nGrid; n++){
      VectorD3 grid = _startPoint + delta * n * _normalizedDirection;
      _gridList[n] = grid;
    }    
  }
  
  void setConnectivityCellGrids(const Array<VectorD3>& cellCentroid)
  {
    /*** method-1 find the closest grid point to each cell
	 so each cell is connected to one grid
	 but does not gurantee each grid has one or more cells
	 if a grid point does not connect to any cell
	 then its conduction_band will be zero due to no interpolation
    ***/
    
    _cellGrids = shared_ptr<CRConnectivity>( new CRConnectivity( _cells, _grids));
    (*_cellGrids).initCount();
    
    for(int c=0; c<_cellList.getLength(); c++){
      const int cellIndex = _cellList[c];
      (*_cellGrids).addCount(cellIndex, 1);
    }
    (*_cellGrids).finishCount();
    
    for(int c=0; c<_cellList.getLength(); c++){
      const int cellIndex = _cellList[c];
      VectorD3 proj = projectionFromPointToLine(_startPoint, _endPoint, cellCentroid[cellIndex]);
      double shortest = 1.0e40;
      int gridIndex = 0;
      for(int g=0; g<_grids.getCount(); g++){
	double distance = mag(proj - _gridList[g]);
	if (distance < shortest){
	  shortest = distance;
	  gridIndex = g;
	}
      }    
      (*_cellGrids).add(cellIndex, gridIndex);
    }
    (*_cellGrids).finishAdd();

 

    /*** method-2 for each grid point, find all the cells within a certain range
	 as long as the range is big enough, the connection is ok
	 but the search range has to be setup mannually
    ***/
    _gridCells = shared_ptr<CRConnectivity>( new CRConnectivity( _grids, _cells));

    const T searchRadius = 0.1e-6/20.0;
    //const T searchRadius = 0.06;

    (*_gridCells).initCount();
    
    const int nGrids = _grids.getCount();
    
    for (int g=0; g < nGrids; g++){
      int nb = 0;
      for (int c=0; c < _cellList.getLength(); c++){
	const int cellIndex = _cellList[c];
	T distance = mag( _gridList[g] - cellCentroid[cellIndex]);
	if (distance < searchRadius){
	  nb++;
	}
      }
      (*_gridCells).addCount(g, nb);
    }
    
    (*_gridCells).finishCount();
    
    for (int g=0; g < nGrids; g++){
      for (int c=0; c < _cellList.getLength(); c++){
	const int cellIndex = _cellList[c];
	T distance = mag( _gridList[g] - cellCentroid[cellIndex]);
	if (distance < searchRadius){
	  (*_gridCells).add(g, cellIndex);
	}
      } 
    }

    (*_gridCells).finishAdd(); 
  
  }

  // interpolation matrix from cell to grid
  void calculateInterpolationMatrix(const int method, 
				    const Array<VectorD3>& cellCentroid)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    
    const CRConnectivity& gridCells = *_gridCells; 
    const CRConnectivity& cellGrids = *_cellGrids; 
    
    const Array<int>& gcRow = gridCells.getRow();
    const Array<int>& gcCol = gridCells.getCol();
    
    const Array<int>& cgRow = cellGrids.getRow();
    const Array<int>& cgCol = cellGrids.getCol();
    
    _cellToGrid = shared_ptr<IMatrix> (new IMatrix(gridCells));
    _gridToCell = shared_ptr<IMatrix> (new IMatrix(cellGrids));

    Array<T>& cellToGridCoeff = _cellToGrid->getCoeff();  
    Array<T>& gridToCellCoeff = _gridToCell->getCoeff();  
   
    if(method == 1){
      /*** algebra average interpolation ***/
      const int nGrids = _gridList.getLength(); 
      for(int n = 0; n < nGrids; n++){
	for(int nc=gcRow[n]; nc<gcRow[n+1]; nc++){
	  const int nCells = gcRow[n+1] - gcRow[n];
	  T wt = 1.0 / nCells;
	  cellToGridCoeff[nc] = wt;
	}
      }
      /*
      const int nLocalCells = _cellList.getLength();
      for(int c=0; c < nLocalCells; c++){
	const int cellIndex = _cellList[c];
	for (int nc = cgRow[cellIndex]; nc < cgRow[cellIndex+1]; nc++){
	  const int nGrids = cgRow[cellIndex+1] - cgRow[cellIndex];
	  T wt = 1.0 / nGrids;
	  gridToCellCoeff[nc] = wt;
	}
      }
      */
    }

    if(method == 2){
      /*** distance weighted interpolation ***/
      const int nGrids = _gridList.getLength(); 
      for(int n = 0; n < nGrids; n++){
	T wtSum(0);
	int nnb(0);
	for(int nc=gcRow[n]; nc<gcRow[n+1]; nc++){
	  const int c = gcCol[nc];
	  VectorT3 dr(cellCentroid[c] - _gridList[n]);
	  T wt = 0;
	  if (dot(dr,dr) < 1e-50) {
	    wt = 1.0e10;
	  }
	  else {
	    wt = 1.0/dot(dr,dr);
	  }
	  cellToGridCoeff[nc] = wt;
	  wtSum += wt;
	  nnb++;
	}
	
	if (nnb == 0)
	  throw CException("no cell neighbors for grid");
	
	for(int nc=gcRow[n]; nc<gcRow[n+1]; nc++){
	  cellToGridCoeff[nc] /= wtSum;
	}
      }
      /*
      const int nLocalCells = _cellList.getLength();
      for(int c=0; c < nLocalCells; c++){
	T wtSum(0);
	T nnb(0);
	const int cellIndex = _cellList[c];
	for(int nc = cgRow[cellIndex]; nc < cgRow[cellIndex+1]; nc++){
	  const int g = cgCol[nc];
	  VectorT3 dr(cellCentroid[cellIndex] - _gridList[g]);
	  T wt = 0;
	  if (dot(dr,dr) < 1e-50) {
	    wt = 1.0e10;
	  }
	  else {
	    wt = 1.0/dot(dr,dr);
	  }
	  gridToCellCoeff[nc] = wt;
	  wtSum += wt;
	  nnb++;
	}

	//if (nnb == 0)
	//throw CException("no grid neighbors for cell");
	
	for(int nc = cgRow[cellIndex]; nc < cgRow[cellIndex+1]; nc++){
	  gridToCellCoeff[nc] /= wtSum;
	}
      }
      */
    }

    if(method == 3){
      /*** linear least square interpolation ***/
    }
  }


  
  void interpolateConductionBand(const Array<T>&  cellValue)
  {
    _conduction_band.zero();
    
    _cellToGrid->multiplyAndAdd(_conduction_band, cellValue);
    
   }

  
  void interpolateTransmission(Array<T>&  cellValue)
  {
    const int nLocalCells = _cellList.getLength();
    for(int c=0; c<nLocalCells; c++){
      const int cellIndex = _cellList[c];
      for(int g=0; g<_cellGrids->getCount(cellIndex); g++){
	const int gridIndex = (*_cellGrids)(cellIndex, g);
	cellValue[cellIndex] =  _transmission[gridIndex];
      }
    }
    
    //cellValue.zero();
    //_gridToCell->multiplyAndAdd(cellValue, _transmission);
  
  }

  
  void calculateTransmission(const T& energy,
			     const T& dielectric_ionization,
			     const T& electron_effmass,
			     const string& flag)
  {
    const int nGrids = _gridList.getLength();
   
    _transmission.zero();

    if (flag == "substrate"){

      _transmission[0] = 1.0;

      const T factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
      //const T factor = -1.0;
      for(int g=1; g<nGrids; g++){
	T delta = mag( _gridList[g] - _gridList[g-1]);
	T valueI = PositiveValueOf( _conduction_band[g] - energy);
	T valueJ = PositiveValueOf( _conduction_band[g-1] - energy);
	T avg = (valueI + valueJ) / 2.0;
	T exponent = factor * sqrt(avg) * delta;
	_transmission[g] = _transmission[g-1] * exp(exponent);
      }   
    }

    else if (flag == "membrane"){
     
      _transmission[nGrids-1] = 1.0;
      
      const T factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
      //const T factor = -1.0;
      for(int g=nGrids-2; g>=0; g--){
	T delta = mag( _gridList[g+1] - _gridList[g]);
	T valueI = PositiveValueOf( _conduction_band[g+1] - energy);
	T valueJ = PositiveValueOf( _conduction_band[g] - energy);
	T avg = (valueI + valueJ) / 2.0;
	T exponent = factor * sqrt(avg) * delta;
	_transmission[g] = _transmission[g+1] * exp(exponent);;
      }   
    }

    else 
      throw CException ("unknown transmission calculation flag!");
  }

 
  void outputGrids(string file)
  {
    char *f = &file[0];
    FILE *fp = fopen(f, "a");
    const int nGrid = _grids.getCount();
    for(int n=0; n<nGrid; n++){
      fprintf(fp, "%i\t%e\t%e\t%e\n", n, _gridList[n][0], _gridList[n][1], _gridList[n][2]);
    }
    fclose(fp);
  }

  void outputCenterLine(string file)
  { 
    char *f = &file[0];
    FILE *fp = fopen(f, "a");
    fprintf(fp, "start %e\t%e\t%e\n", _startPoint[0],_startPoint[1],_startPoint[2]); 
    fprintf(fp, "end   %e\t%e\t%e\n", _endPoint[0],_endPoint[1], _endPoint[2]);
    fprintf(fp, "direction   %e\t%e\t%e\n", _normalizedDirection[0],_normalizedDirection[1], _normalizedDirection[2]);
  }
  
  void outputCellList(string file)
  {
    char *f = &file[0];
    FILE *fp = fopen(f, "a"); 
    fprintf(fp, "I own the following cells\n");
    for(int c=0; c<_cellList.getLength(); c++){
      fprintf(fp, "%i\t%i\n", c, _cellList[c]); 
    }
    fclose(fp);
  }
  
  void outputTransmission(string file)
  {
    char *f = &file[0];
    FILE *fp = fopen(f, "a"); 
    fprintf(fp, "here is the transmission on grids\n");
    const int nGrid = _grids.getCount();
    for(int n=0; n<nGrid; n++){
      fprintf(fp, "%i\t%e\t%e\n", n, _conduction_band[n], _transmission[n]);
    }
    fclose(fp);
  }

  void outputConnectivityCellGrids(string file)
  {
    char *f = &file[0];
    FILE *fp = fopen(f, "a");
    fprintf(fp, "the following is connectivity from cell to grid\n");
    for(int c=0; c<_cellList.getLength(); c++){
      const int cellIndex = _cellList[c];
      for(int g=0; g<_cellGrids->getCount(cellIndex); g++){
	const int gridIndex = (*_cellGrids)(cellIndex, g);
	fprintf(fp, "%i\t%i\n", cellIndex, gridIndex);
      }
    }

    fprintf(fp, "the following is connectivity from grid to cell\n");
    for(int g=0; g<_gridList.getLength(); g++){
      for(int c=0; c<_gridCells->getCount(g); c++){
	const int cellIndex = (*_gridCells)(g, c);
	fprintf(fp, "%i\t%i\n", g, cellIndex);
      }
    }

    fclose(fp);
  }
 


 private:
  DielectricOneDimColumn<T> (const DielectricOneDimColumn<T>&);
  VectorD3 _startPoint;
  VectorD3 _endPoint;
  VectorD3 _normalizedDirection;
  T _radius;
  StorageSite _grids;
  const StorageSite& _cells;

  Array<VectorD3> _gridList;
  Array<int> _cellList;
  Array<T>  _conduction_band;
  Array<T>  _transmission;    

  shared_ptr<CRConnectivity> _cellGrids;
  shared_ptr<CRConnectivity> _gridCells;
  shared_ptr <CRMatrixTranspose<T,T,T> > _cellToGrid;
  shared_ptr <CRMatrixTranspose<T,T,T> > _gridToCell;

 

};

typedef Vector<double, 3> VecD3;

template<class T> 
void setCenterLines(const int nXCol, const int nYCol, 
		       const VecD3 corner1_1, const VecD3 corner1_2, 
		       const VecD3 corner1_3, const VecD3 corner1_4,
		       const VecD3 corner2_1, const VecD3 corner2_2, 
		       const VecD3 corner2_3, const VecD3 corner2_4,
		       vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList, 
		       StorageSite& columns)
{
  typedef Vector<double,3> VectorD3;
  

  double xFrom, xTo, yFrom, yTo, zFrom, zTo, deltaX, deltaY, size=0;
  //mesh the boundary surface 
  /* here we assume x-y plane is the boundary surface
     z direction is the electron transmission direction 
     four corner points define the geometry of the surface 
  */
    
  xFrom = findMin(corner1_1[0],corner1_2[0],corner1_3[0],corner1_4[0]);
  xTo   = findMax(corner1_1[0],corner1_2[0],corner1_3[0],corner1_4[0]);
  yFrom = findMin(corner1_1[1],corner1_2[1],corner1_3[1],corner1_4[1]);
  yTo   = findMax(corner1_1[1],corner1_2[1],corner1_3[1],corner1_4[1]);
  zFrom = findMin(corner1_1[2],corner1_2[2],corner1_3[2],corner1_4[2]);
  zTo   = findMax(corner1_1[2],corner1_2[2],corner1_3[2],corner1_4[2]);

  deltaX = (xTo - xFrom) / nXCol;
  deltaY = (yTo - yFrom) / nYCol;

  if(deltaX > size) size = deltaX;
  if(deltaY > size) size = deltaY;
 
  for ( int i=0; i<nXCol; i++){
    for ( int j=0; j<nYCol; j++){
      VectorD3 point; 
      point[0] = xFrom + deltaX * (i+1) - deltaX/2.0;
      point[1] = yFrom + deltaY * (j+1) - deltaY/2.0;
      point[2] = (zTo + zFrom)/2.0;
      int col = i*nXCol + j;
      columnList[col]->setStartPoint(point);
    }
  }
  
  xFrom = findMin(corner2_1[0],corner2_2[0],corner2_3[0],corner2_4[0]);
  xTo   = findMax(corner2_1[0],corner2_2[0],corner2_3[0],corner2_4[0]);
  yFrom = findMin(corner2_1[1],corner2_2[1],corner2_3[1],corner2_4[1]);
  yTo   = findMax(corner2_1[1],corner2_2[1],corner2_3[1],corner2_4[1]);
  zFrom = findMin(corner2_1[2],corner2_2[2],corner2_3[2],corner2_4[2]);
  zTo   = findMax(corner2_1[2],corner2_2[2],corner2_3[2],corner2_4[2]);
  
  deltaX = (xTo - xFrom) / nXCol;
  deltaY = (yTo - yFrom) / nYCol;

  if(deltaX > size) size = deltaX;
  if(deltaY > size) size = deltaY;

  for ( int i=0; i<nXCol; i++){
    for ( int j=0; j<nYCol; j++){
      VectorD3 point; 
      point[0] = xFrom + deltaX * (i+1) - deltaX/2.0;
      point[1] = yFrom + deltaY * (j+1) - deltaY/2.0;
      point[2] = (zTo + zFrom)/2.0;
      int col = i*nXCol + j;
      columnList[col]->setEndPoint(point);
    }
  }  
  
  //set up the normalized direction for all center lines
  const int nColumn = columnList.size();
  for(int n=0; n<nColumn; n++){
    const VectorD3 start = columnList[n]->getStartPoint();
    const VectorD3 end = columnList[n]->getEndPoint();
    VectorD3 nd = (end - start) / (mag(end-start));
    columnList[n]->setNormalizedDirection(nd);
    columnList[n]->setRadius(size);
  }
}

/* set up the connectivity from cell to columns
   assume each cell only connects to one column
   to which it has the shortest distance
*/

template<class T>
void
setConnectivityCellColumns(const StorageSite& cells,
			   const StorageSite& columns,  
			   const Array<Vector<T,3> >& cellCentroid, 
			   vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList, 	  
			   shared_ptr<CRConnectivity>& cellColumns)
{
  const int nColumns = columnList.size();
  
  const int nCells = cells.getCount();
  
  cellColumns = shared_ptr<CRConnectivity>(new CRConnectivity (cells, columns));

  (*cellColumns).initCount();       //step 1//

  for(int c=0; c<nCells; c++){
    (*cellColumns).addCount(c, 1);   //step 2//
  }
	
  (*cellColumns).finishCount();     //step 3//
	
  int columnIndex = 0;
  for(int c=0; c<nCells; c++){
    T shortest = 10000000000;
    for (int n=0; n<nColumns; n++){
      const VectorD3 start = columnList[n]->getStartPoint();
      const VectorD3 end = columnList[n]->getEndPoint();
      T distance = distanceFromPointToLine(start, end, cellCentroid[c]);
      if (distance < shortest){
	shortest = distance;
	columnIndex = n;
      }
    }
    
    (*cellColumns).add(c, columnIndex);   //step 4//
  }
	      
  (*cellColumns).finishAdd();     //step 5//
  

}


void outputConnectivityCellColumns(const Mesh& mesh, const Array<VectorD3>& cellCentroid,
				   shared_ptr<CRConnectivity>& cellColumns)
{
  //check if the connectivity is right
  string fileName = "/home/lin/work/app-memosa/src/fvm/verification/ElectricModel/tunneling/cellToColumns.dat";
  char * file;
  file = &fileName[0];
  FILE *fp = fopen(file, "w");

  const StorageSite& cells = mesh.getCells();
  const int nCells = cells.getCount();
  for (int c=0; c<nCells; c++){
    const int size = (*cellColumns).getCount(c);
    for ( int n=0; n<size; n++){
      int col = (*cellColumns)(c, n);
      fprintf(fp, "%i\t %e\t %e\t %e\t %i\n", c, cellCentroid[c][0],cellCentroid[c][1],cellCentroid[c][2],col); 
    }
  }
  fclose(fp);
  
}



template<class T>
void discretizeCenterLine(vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList)
{
  const int nCol = columnList.size();
  for(int n=0; n<nCol; n++){
    columnList[n]->discretizeCenterLine();
  }
}


template <class T>
void setCellList(vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList, 
		 shared_ptr<CRConnectivity>& columnCells)
{
  const int nCol = columnList.size();
  for(int col=0; col<nCol; col++){
    const int nCells = columnCells->getCount(col);
    Array<int>& cellList = columnList[col]->getCellList();
    cellList.resize(nCells);
    for(int c=0; c<nCells; c++){
      const int cellIndex = (*columnCells)(col, c);
      columnList[col]->setCellList(c, cellIndex);
    }
  }
}


template <class T>
void setConnectivityCellGrids(vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList,
			      const Array<VectorD3>& cellCentroid)
{
  const int nCol = columnList.size();
  for(int col=0; col<nCol; col++){
    columnList[col]->setConnectivityCellGrids(cellCentroid);
  }
}
  
template<class T>
void outputColumn( vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList)
{
  string fileName = "/home/lin/work/app-memosa/src/fvm/verification/ElectricModel/tunneling/columns.dat";
  const int nCol = columnList.size();
  for(int n=0; n<1; n++){
     columnList[n]->outputCenterLine(fileName);
     columnList[n]->outputGrids(fileName);
     columnList[n]->outputCellList(fileName);
    columnList[n]->outputConnectivityCellGrids(fileName);
    columnList[n]->outputTransmission(fileName);
  }
}  


template<class T>
void calculateInterpolationMatrix(vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList,
				  const int method, 
				  const Array<VectorD3>& cellCentroid)
{
   const int nCol = columnList.size();
   for(int n=0; n<nCol; n++){
     columnList[n]->calculateInterpolationMatrix(method, cellCentroid);
   }
}


template<class T>
void ElectronTransmissionCoefficient 
              (const T& energy, 
	       Array<T>& transmission,
	       const Array<T> & conduction_band,
	       const T& dielectric_ionization,
	       const T& electron_effmass,
	       const string& flag,
	       vector<shared_ptr<DielectricOneDimColumn<T> > >& columnList)
{
  const int nCol = columnList.size();

  for(int n=0; n<nCol; n++){

    //interpolate conduction band from cells to columns
    columnList[n]->interpolateConductionBand(conduction_band);
   
    //each column calculate transmission and then interpolate back to cells
    columnList[n]->calculateTransmission(energy, dielectric_ionization, electron_effmass, flag );
   
    columnList[n]->interpolateTransmission(transmission);
  }
}

#endif
