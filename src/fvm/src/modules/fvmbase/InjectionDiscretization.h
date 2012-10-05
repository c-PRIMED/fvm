// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _INJECTIONDISCRETIZATION_H_
#define _INJECTIONDISCRETIZATION_H_

#include "PhysicsConstant.h"
#include "ElectricBC.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "CRMatrix.h"

/**************************************************
Diag type: 2x2Tensor
      | d00,  d01 |
      | d10,  d11 |    

OffDiag type: 2x2Tensor
      | o00,  o01 |
      | o10,  o11 |    

X type: VectorT2
      | x0 |
      | x1 |         

"0" is trapped charge
"1" is band charge

Injection model only modifies x1 

*************************************************/

template <class X, class Diag, class OffDiag>
class InjectionDiscretization : public Discretization
{

 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;
  typedef Array<Vector<T_Scalar, 3> > VectorT3Array;


  InjectionDiscretization(const MeshList& meshes,
			  const GeomFields& geomFields,
			  const Field& varField,
			  const Field& electricField,
			  const Field& conductionbandField,
			  const ElectricModelConstants<T_Scalar>& constants
			  ) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _electricField(electricField),
    _conductionbandField(conductionbandField),
    _constants(constants)
     {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const int nCells = cells.getSelfCount();

    //const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array&> (_electricField[cells]);
    
    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    TArray* ts = (new TArray(cells.getCount()));
    *ts = 0;
    TArray& transmission = *ts;

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    //CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
     
    const VectorT3Array& cellCentroid = 
      dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);
    
    //const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    //DiagArray& diag = matrix.getDiag();

    //OffDiagArray& offdiag = matrix.getOffDiag();

    const T_Scalar dielectric_thickness = _constants["dielectric_thickness"];
    const T_Scalar electron_effmass = _constants["electron_effmass"];
    const T_Scalar temperature = _constants["OP_temperature"];
    //const T_Scalar substrate_voltage = _constants["substrate_voltage"];
    //const T_Scalar membrane_voltage = _constants["membrane_voltage"];
    const T_Scalar fermilevelsubstrate = -_constants["substrate_workfunction"] - _constants["substrate_voltage"];
    //const T_Scalar fermilevelmembrane = -_constants["substrate_workfunction"] - _constants["membrane_voltage"];
    const int subID = _constants["substrate_id"];
    //const int memID = _constants["membrane_id"];
    const int nLevel = _constants["nLevel"];
    const int normal = _constants["normal_direction"];
    const int nTrap = _constants["nTrap"];

    T_Scalar fluxCoeff(0), fermilevel(0); // scatterfactor(0);
    //T_Scalar sourceInjection(0);
    const T_Scalar alpha = 4.0 * PI * (electron_effmass*ME) / pow(H_SI, 3.0);
    //const T_Scalar energystep =  fabs(substrate_voltage - membrane_voltage) / nLevel;
    const T_Scalar energystep = 0.01;
    //=======================================//
    // injection from substrate to dielectric 
    //=======================================//
    for(int c=0; c<cells.getCount(); c++){
      transmission[c] = 0.0;
    } 

    fermilevel = fermilevelsubstrate;
   
    for (T_Scalar en = fermilevel-4.0; en <= fermilevel+4.0; en += energystep){   
   
      const T_Scalar supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);

      const T_Scalar fermifunction = FermiFunction(en, fermilevel, temperature);

      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;

	if (fg.id == subID){
	  const StorageSite& faces = fg.site;	 
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const CRConnectivity& cellCells = mesh.getCellCells();
	  const int nFaces = faces.getCount();
	  
	  for(int f=0; f<nFaces; f++){

	    //--------------------------------------------------------------------------------------------
	    // 1D int array to store the cell index along each straight line from substrate to membrane
	    int indices[nLevel+1];    

	    int c0 = faceCells(f,0);
	    int c1 = faceCells(f,1);

	    indices[0] = c1;
	    	     
	    int low = c1;
	    int me = c0;	      
	    int high = c0;
	       
	    for(int l=0; l < nLevel; l++){
	      if (me > nCells) 
		throw CException("index out of boundary in elec_model injection calculation");
	      indices[l+1] = me;
	      const int nbc = cellCells.getCount(me);
	      T_Scalar drmin = 0.0;
	      int neighborUp = 0;
	      for(int nc = 0; nc < nbc; nc++){
		const int neighbor = cellCells(me, nc);
		const T_Scalar dr = cellCentroid[me][normal] - cellCentroid[neighbor][normal];
		if (dr < drmin){
		  drmin = dr;
		  neighborUp = neighbor;
		}
	      }	      
	      high = neighborUp;
	      low = me;
	      me = high;	 
	    }
	    //--------------------------------------------------------------------------------------------
	    //calculate transmission coefficient along the line
	    transmission[indices[0]]=1.0;
	    for(int l=0; l < nLevel; l++){
	      int low = indices[l];
	      int me = indices[l+1];
	      
	      T_Scalar dX = cellCentroid[me][normal] - cellCentroid[low][normal];
	      T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	      T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
	      T_Scalar valueLow = PositiveValueOf( conduction_band[low] - en);
	      T_Scalar avg = (valueMe + valueLow) / 2.0;
	      T_Scalar exponent = factor * sqrt(avg) * fabs(dX);
#if DEBUG
	      (*mark)[me] = 1;
#endif		
	      transmission[me] = transmission[low] * exp(exponent);
	    }
	    //--------------------------------------------------------------------------------------------
	    for(int l=0; l < nLevel; l++){
	      //int low = indices[l];
	      int me = indices[l+1];
	      //injection occurs where the sign of (en-conductionband) switches. 
	      if ((en-conduction_band[me]) > 0) {
		
		const T_Scalar dX = dielectric_thickness/nLevel;		
		fluxCoeff = alpha * transmission[me] * supplyfunction * fermifunction * energystep * QE * cellVolume[me] / fabs(dX) ;  
		rCell[me][nTrap] += fluxCoeff; 
		break;
	      }
	     
	    }
	  }
	}
      }
    }


#if 0
    //=======================================//
    // Injection from membrane to dielectric 
    //=======================================//
    for(int c=0; c < cells.getCount(); c++){
      transmission[c] = 0.0;
    }

    fermilevel = fermilevelmembrane;
    for (T_Scalar en = fermilevel-4.0; en <= fermilevel+4.0; en += energystep){   
   
      const T_Scalar supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);

      const T_Scalar fermifunction = FermiFunction(en, fermilevel, temperature);

      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;

	if (fg.id == memID){
	  const StorageSite& faces = fg.site;	 
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const CRConnectivity& cellCells = mesh.getCellCells();
	  const int nFaces = faces.getCount();
	  
	  for(int f=0; f<nFaces; f++){

	    //--------------------------------------------------------------------------------------------
	    // 1D int array to store the cell index along each straight line from substrate to membrane
	    int indices[nLevel+1];    

	    int c0 = faceCells(f,0);
	    int c1 = faceCells(f,1);

	    indices[0] = c1;
	    	     
	    int low = c1;
	    int me = c0;	      
	    int high = c0;
	       
	    for(int l=0; l < nLevel; l++){
	      if (me > nCells) 
		throw CException("index out of boundary in elec_model injection calculation");
	      indices[l+1] = me;
	      const int nbc = cellCells.getCount(me);
	      T_Scalar drmin = 0.0;
	      int neighborUp = 0;
	      for(int nc = 0; nc < nbc; nc++){
		const int neighbor = cellCells(me, nc);
		const T_Scalar dr = cellCentroid[me][normal] - cellCentroid[neighbor][normal];
		if (dr > drmin){
		  drmin = dr;
		  neighborUp = neighbor;
		}
	      }	      
	      high = neighborUp;
	      low = me;
	      me = high;	 
	    }
	    //--------------------------------------------------------------------------------------------
	    //calculate transmission coefficient along the line
	    transmission[indices[0]]=1.0;
	    for(int l=0; l < nLevel; l++){
	      int low = indices[l];
	      int me = indices[l+1];
	      
	      T_Scalar dX = cellCentroid[me][normal] - cellCentroid[low][normal];
	      T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	      T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
	      T_Scalar valueLow = PositiveValueOf( conduction_band[low] - en);
	      T_Scalar avg = (valueMe + valueLow) / 2.0;
	      T_Scalar exponent = factor * sqrt(avg) * fabs(dX);
#if DEBUG
	      (*mark)[me] = 1;
#endif		
	      transmission[me] = transmission[low] * exp(exponent);
	    }
	    //--------------------------------------------------------------------------------------------
	    for(int l=0; l < nLevel; l++){
	      int low = indices[l];
	      int me = indices[l+1];
	      //injection occurs where the sign of (en-conductionband) switches. 
	      if ((en-conduction_band[me]) > 0) {
		
		const T_Scalar dX = energystep / electric_field[me][normal];
	      
		fluxCoeff = alpha * transmission[me] * supplyfunction * fermifunction * energystep * QE * cellVolume[me] / fabs(dX) ;
	  
		rCell[me][nTrap] += fluxCoeff; 
		//cout << "transmission " << me << " " << transmission[me] << endl;
		
		break;
	      }
	    }
	  }
	}
      }
    }
    
#endif


  }
 
	


 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _electricField;
  const Field& _conductionbandField;
  const ElectricModelConstants<T_Scalar>& _constants;
 

};





#endif
