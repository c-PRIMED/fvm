// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CONVECTIONDISCRETIZATION_KMODEL_H_
#define _CONVECTIONDISCRETIZATION_KMODEL_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include <math.h>
#include "Discretization.h"
#include "StorageSite.h"
#include "Gradient.h"
#include <cassert>

template<class X, class Diag, class OffDiag>
class ConvectionDiscretization_Kmodel : public Discretization
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Gradient<X> XGrad;

  typedef Array<int> IntArray;
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;

  ConvectionDiscretization_Kmodel(const MeshList& meshes,
				  const GeomFields& geomFields,
				  Field& varField,
				  const double cx,
				  const double cy,
				  const double cz,
				  bool useCentralDifference=false,
                                  bool useMinmodFluxLimiter=true) :
    //const bool useCentralDifference=false) :
  Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField), 
    _cx(cx),
    _cy(cy),
    _cz(cz),
    _useCentralDifference(useCentralDifference),
    _useMinmodFluxLimiter(useMinmodFluxLimiter),
    _cellCells ( (meshes[0])->getCellCells() ),
    _cellCells2( (meshes[0])->getCellCells2() ),
    _faceCells ( (meshes[0])->getAllFaceCells() )
  {}
    
 
  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {

    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();


    // should there be some other checks ?
    //if (!_convectingFluxField.hasArray(faces))
    //  return;

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix =
      dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));

    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const CRConnectivity& cellCells = mesh.getCellCells();
    const CRConnectivity& cellCells2 = mesh.getCellCells2();    

    CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
    DiagArray& diag = matrix.getDiag();

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    //const GradArray& xGradCell = dynamic_cast<GradArray>(_varGradientField[cells]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
    
    const int nFaces = faces.getCount();

    const VectorT3Array& faceArea  = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);   
    const VectorT3Array& cellCoord = dynamic_cast<const VectorT3Array& >(_geomFields.coordinate[cells]); 
    _cellCoord = &cellCoord;
    //const X nondim_length=_options["nonDimLt"];
    //const X Lx=_options["nonDimLx"];
    //const X Ly=_options["nonDimLy"];
    //const X Lz=_options["nonDimLz"];
    
    if (_geomFields.gridFlux.hasArray(faces)){
        shared_ptr<TArray> gridFluxPtr(new TArray(nFaces));
	TArray& gridFlux = *gridFluxPtr;
        gridFlux = dynamic_cast<const TArray&>(_geomFields.gridFlux[faces]);

	for(int f=0; f<nFaces; f++){
            const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);
	    //const T_Scalar faceCFlux = convectingFlux[f] - gridFlux[f];
	    const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz - gridFlux[f]; 
	    //const T_Scalar faceCFlux = faceArea[f][0]*_cx*nondim_length/Lx+faceArea[f][1]*_cy*nondim_length/Ly+faceArea[f][2]*_cz*nondim_length/Lz - gridFlux[f];

	    X varFlux;
	    if (faceCFlux > T_Scalar(0)){
	        varFlux = faceCFlux*xCell[c0];
                diag[c0] -= faceCFlux;
                assembler.getCoeff10(f) += faceCFlux;
	    }
	    else{
                varFlux = faceCFlux*xCell[c1];
                diag[c1] += faceCFlux;
                assembler.getCoeff01(f)-= faceCFlux;
	    }
        
	    rCell[c0] -= varFlux;
	    rCell[c1] += varFlux;
	}
    }  else {
      if (_useCentralDifference){
	for(int f=0; f<nFaces; f++){
              const int c0 = faceCells(f,0);
              const int c1 = faceCells(f,1);
	      const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz;
              //const T_Scalar faceCFlux = faceArea[f][0]*_cx*nondim_length/Lx+faceArea[f][1]*_cy*nondim_length/Ly+faceArea[f][2]*_cz*nondim_length/Lz;
              bool isIBFace = (((ibType[c0] == Mesh::IBTYPE_FLUID)
                                && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
                               ((ibType[c1] == Mesh::IBTYPE_FLUID)
                                && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)));
            

              X varFlux =0.5*faceCFlux*(xCell[c0] + xCell[c1]);

              rCell[c0] -= varFlux;
              rCell[c1] += varFlux;

              if (isIBFace){
                  // linearize the actual flux as calculated
                  // above. this will ensure that the Ib
                  // discretization will be able to fix the value
                  // correctly using the ib face value
                  
                  diag[c0] -= 0.5*faceCFlux;
                  assembler.getCoeff10(f) -= 0.5*faceCFlux;
                  diag[c1] += 0.5*faceCFlux;
                  assembler.getCoeff01(f) += 0.5*faceCFlux;
              } else {
                  // linearize as upwind flux so that linear system
                  // remains diagonally dominant
                  if (faceCFlux > T_Scalar(0)){
                      diag[c0] -= faceCFlux;
                      assembler.getCoeff10(f) += faceCFlux;
                  } else {
                      diag[c1] += faceCFlux;
                      assembler.getCoeff01(f)-= faceCFlux;
                  }
              }
          }
      } /*else if (_useMinmodFluxLimiter){    ///////////////////////////////////////////Edit by Mizuki (Aug 2nd,2012)
	for (int f=0; f<nFaces; f++){
	  const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);

          const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz;
        //'qsub -I -q  prism -l select=1:ncpus=1,walltime=240:00:00 -v DISPLAY' 
          X varFlux;

	  // Need definition of flux limiter
	  const int flux_east = 1.0;
          const int flux_west = 1.0;
          // CD=cental difference, LUD=linear upwind difference, UP=upwind difference
	  bool CD_CD   = false;
          bool CD_LUD  = false;
	  bool CD_UD   = false;
	  bool LUD_CD  = false;
	  bool LUD_LUD = false;
	  bool LUD_UD  = false;
  	  bool UD_CD   = false;
	  bool UD_LUD  = false;
  	  bool UD_UD   = false;
	  
	  if((flux_east >= 1.0) && (flux_west >= 1.0)){ 
	      bool CD_CD = true; // case1   
	  }
	  else if ((flux_east >= 1.0) && (flux_west > 0.0) && (flux_west < 1.0)){
              bool CD_LUD = true; // case2
	  }
	  else if ((flux_east >= 1.0) && (flux_west <= 0.0)){
	      bool CD_UD = true;  // case 3
	  }
	  else if ((flux_east > 0.0) && (flux_east < 1.0) && (flux_west >= 1.0)){
              bool LUD_CD = true; // case4
	  }
	  else if ((flux_east > 0.0) && (flux_east < 1.0) && (flux_west > 0.0) && (flux_west< 1.0)){
	      bool LUD_LUD = true; // case5
	  }
	  else if ((flux_east > 0.0) && (flux_east < 1.0) && (flux_west <= 0)){
	      bool LUD_CD = true; // case6
          }
	  else if ((flux_east <= 0.0) && (flux_west >= 1.0)){
	      bool UD_CD = true; // case7
	  }
	  else if ((flux_east <= 0.0) && (flux_west > 0.0) && (flux_west < 1.0)){
	      bool UD_LUD = true; // case8
	  }
	  else if ((flux_east <= 0.0) && (flux_west <= 0.0)){
	      bool UD_UD = true; // case9
	  }

	  
	  if(CD_CD){ // case1
	      if (faceCFlux > T_Scalar(0)){
	          const int c11 = identifyNeighborCell(f,  2, faceCFlux, true);
		  diag[c0] -= 0.5*faceCFlux;  // (i-1)
	          assembler.getCoeff10(f) += 0; // (i)
	  	  assembler.getCoeff02(f) += 0.5*faceCflux; // (i+1)

		  rCell[c0] -= 0.5*faceCFlux*xCell[c0];
		  rCell[c1] += 0;
		  rCell[c11] += 0.5*faceCFlux*xCell[c11]; 
		  }
	      else{
	          const int c00 = identifyNeighborCell(f, 2, faceCFlux, true); 
		  diag[c1] += 0.5*faceCFlux; // (i+1)
		  assembler.getCoeff01(f) += 0; // (i)
		  assembler.getCoeff20(f) -= 0.5*faceCFlux; // (i-1)

		  rCell[c0] -= 0;
		  rCell[c1] += 0.5*faceCFlux*xCell[c1];
		  rCell[c00] -= 0.5*faceCFlux*xCell[c00];	
		  }
	       }
	  else if (CD_LUD){ // case2 
	      if(faceCFlux > T_Scalar(0)){
                  c00 = identifyNeighborCell(f, 1, faceCFlux, false);
	          c11 = identifyNeighborCell(f, 2, faceCFlux, true);    
	          diag[c0] -= 1.5*faceCFlux; // (i-1)
		  assembler.getCoeff10(f) += 0.5*faceCFlux; // (i)
	          assembler.getCoeff220(f) += 0.5*faceCFlux; // (i-2)
	          assembler.getCoeff20(f) += 0.5*faceCFlux; // (i+1)
	
		  rCell[c0] -= 1.5*faceCFlux*xCell[c0];
		  rCell[c1] += 0.5*faceCFlux*xCell[c1];
		  rCell[c00] += 0.5*faceCFlux*xCell[c00];
	 	  rCell[c11] += 0.5*faceCFlux*xCell[c11];
	       }
	      else{
		  diag[c1] += faceCFLux;  // (i+1)
		  assembler.getCoeff01(f) -= faceCFlux; // (i)
	          
	          rCell[c0] -= faceCFlux*xCell[c0];
		  rCell[c1] += faceCFlux*xCell[c1];
	      } 
	  }
	  else if (CD_UD){ // case3
              if(faceCFlux > T_Scalar(0)){
		  c11 = identifyNeighborCell(f, 2, faceCFlux, true);
	          diag[c0] -= faceCFlux; // (i-1)
		  assembler.getCoeff10(f) += 0.5*faceCFlux; //(i)
		  assembler.getCoeff20(f) += 0.5*faceCFlux; //(i+1)

		  rCell[c0] -= faceCFlux*xCell[c0];
		  rCell[c1] += 0.5*faceCFluc*xCell[c1];
		  rCell[c11] += 0.5*faceCFlux*xCell[c11];     
	      }
	      else{
	          diag[c1] += 0.5*faceCFlux; // (i+1)
		  assembler.getCoeff01(f) -= 0.5*faceCFlux; //(i)
		  
		  rCell[c0] -= 0.5*faceCFlux*xCell[c0];
		  rCell[c1] += 0.5*faceCFLux*xCell[c1]; 	
	      }		
	  }
	  else if (LUD_CD){ // case 4
	      if(faceCFlux > T_Scalar(0)){
	          diag[c0] -= faceCFlux; // (i-1)
		  assembler.getCoeff10(f) += faceCFlux; // (i)

		  rCell[c0] -= faceCFlux*xCell[c0];
		  rCell[c1] += faceCflux*xCell[c1];
	      }
	      else{
		  c00 = identifyNeighborCell(f, 2, faceCFlux, true);
		  c11 = idenrifyNeighborCell(f, 1, faceCFlux, false);
		  diag[c1] += 1.5*faceCFlux; //(i+1)
		  assembler.getCoeff01(f) -= 0.5*faceCFlux; //(i)
		  assembler.getCoeff02(f) -= 0.5*faceCFlux; // (i-1)
		  assembler.getCoeff022(f) -= 0.5*faceCFlux; // (i+2)

		  rCell[c0] -= 0.5*faceCFlux*xCell[c0];
		  rCell[c1] += 1.5*faceCFlux*xCell[c1];
		  rCell[c00] -= 0.5*faceCFlux*xCell[c00];
		  rCell[c11] -= 0.5*faceCFlux*xCell[c11];   
	      }
	  } 
	  else if (LUD_LUD){ // case 5
              if(faceCFlux > T_Scalar(0)){
	          c00 = identifyNeighborCell(f, 1, faceCFlux, false);
		  diag[c0] -= 2.0*faceCFlux; // (i-1)
		  assembler.getCoeff10(f) += 1.5*faceCFlux; // (i)
		  assembler.getCoeff220(f) += 0.5*faceCFlux; // (i-2)

		  rCell[c0] -= 2.0*faceCFlux*xCell[c0];
		  rCell[c1] += 1.5*faceCFlux*xCell[c1];
		  rCell[c00] += 0.5*faceCFlux*xCell[c00];
	      {
	      else{
	          c11 = identifyNeighborCell(f, 1, faceCFlux, false);
		  diag[c1] += 2*faceCFlux;
		  assembler.getCoeff 
	      }
	  }
	  else if (LUD_UD){i
	  }
	  else if (UD_CD){
	  }
	  else if (UD_LUD){
	  }
	  else if (UD_UD){
	      if(faceCFlux > T_Scalar(0)){
	          varFlux = faceCFlux*xCell[c0];
		  diag[c0] -= faceCFlux;   // (i-1)
		  assembler.getCoeff10(f) += faceCflux; //(i)	
	      }
	      else{
		  varFlux = faceCFlux*xCell[c1];
		  diag[c1] -= faceCFlux;   // (i+1)
		  assembler.getCoeff01(f) -= faceCFlux; // (i)
	      }
	      rCell[c0] -= varFlux;
	      rCell[c1] += varFlux;
	  }


	}
      }*/  else {    //////////////////////////////////////////End of edit by Mizuki
          for(int f=0; f<nFaces; f++) {
              const int c0 = faceCells(f,0);
              const int c1 = faceCells(f,1);
              const T_Scalar faceCFlux = faceArea[f][0]*_cx+faceArea[f][1]*_cy+faceArea[f][2]*_cz;
	      //const T_Scalar faceCFlux = faceArea[f][0]*_cx*nondim_length/Lx+faceArea[f][1]*_cy*nondim_length/Ly+faceArea[f][2]*_cz*nondim_length/Lz;


              X varFlux;
            
              if (faceCFlux > T_Scalar(0)){
                  varFlux = faceCFlux*xCell[c0];
                  diag[c0] -= faceCFlux;
                  assembler.getCoeff10(f) += faceCFlux;
              }  else {
                  varFlux = faceCFlux*xCell[c1];
                  diag[c1] += faceCFlux;
                  assembler.getCoeff01(f)-= faceCFlux;
              }

              rCell[c0] -= varFlux;
              rCell[c1] += varFlux;
           }
      }
  }



    //const int nCells = cells.getSelfCount();
    //for(int c=0;c<nCells;c++)
    //{
    //    const T_Scalar cImb = continuityResidual[c];
    //    diag[c] += cImb;
    //}

  }
private:
  int    identifyNeighborCell(int f, int level, double faceCFlux, const bool direction = true){
	
        double cx = _cx;
        double cy = _cy;
	double cz = _cz;
	int FirstUpwindCell = -1;
	const int c0 = _faceCells(f,0);
	const int c1 = _faceCells(f,1);
	
	
	if (direction == false)
	   faceCFlux = -faceCFlux;
	
	if ( faceCFlux > 0){
	     FirstUpwindCell = c0;
	} else {
	     FirstUpwindCell = c1;  
        }	
	      
	
	 
        if (direction==false){
	   cx = -cx;
	   cy = -cy;
	   cz = -cz;
        }
		 
	
	int c00 = -1;
        double angleCriteria = 2.0*M_PI;
        int ncells2 = _cellCells2.getCount(FirstUpwindCell);
        for (int n=0; n < _cellCells2.getCount(FirstUpwindCell); n++) {
           const int cellID =_cellCells2(FirstUpwindCell,n);
	   if (cellID != c1) {
	      const double vx = (*_cellCoord)[FirstUpwindCell][0] - (*_cellCoord)[cellID][0];
	      const double vy = (*_cellCoord)[FirstUpwindCell][1] - (*_cellCoord)[cellID][1];
	      const double vz = (*_cellCoord)[FirstUpwindCell][2] - (*_cellCoord)[cellID][2];

		    
	      //find angle between ((vx,vy,vz) , c=(cx,cy,cz))
	      double dotProduct = vx*cx + vy*cy + vz*cz;
	      double abs_vc = pow(vx*vx + vy*vy + vz*vz,0.5) * pow(cx*cx + cy*cy + cz*cz, 0.5);
	      double angle = acos(dotProduct/abs_vc);
	      if ( fabs(angle) < angleCriteria){
		 angleCriteria = angle;
		 c00 = cellID;
	      }
	   } 
        }
        assert(c00 != -1); 
        
        return  c00;	    
       
  }
  
  

  const GeomFields& _geomFields;
  const Field& _varField;
  const double _cx;
  const double _cy;
  const double _cz;
  const bool _useCentralDifference;  
  const bool _useMinmodFluxLimiter;
  KineticModelOptions<X> _options;
  const CRConnectivity& _cellCells;
  const CRConnectivity& _cellCells2;
  const CRConnectivity& _faceCells;
  const VectorT3Array* _cellCoord;
};

#endif
