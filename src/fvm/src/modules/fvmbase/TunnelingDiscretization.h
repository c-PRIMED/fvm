// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _TUNNELINGDISCRETIZATION_H_
#define _TUNNELINGDISCRETIZATION_H_

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
#include "ElectricUtilityFunctions.h"


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

Tunneling model only modifies d0 and x0 

*************************************************/


template <class X, class Diag, class OffDiag>
class TunnelingDiscretization : public Discretization
{

 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef Array<int> IntArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;
  typedef Array<Vector<T_Scalar, 3> > VectorT3Array;

  TunnelingDiscretization(const MeshList& meshes,
			  const GeomFields& geomFields,
			  const Field& varField,			  
			  const Field& conductionbandField,
			  const ElectricModelConstants<T_Scalar>& constants, 
			  T_Scalar& fluxIn, 
			  T_Scalar& fluxOut
			  ) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),    
    _conductionbandField(conductionbandField),
      _constants(constants),
      _fluxIn(fluxIn),
      _fluxOut(fluxOut)
      {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    #define DEBUG  0

    const StorageSite& cells = mesh.getCells();

    //const StorageSite& faces = mesh.getFaces();
    
    const int nCells = cells.getSelfCount();      

    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const VectorT3Array& cellCentroid = 
      dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(_varField[cells]);
    
    TArray* ts = new TArray(cells.getCount());
    *ts = 0;
    TArray& transmission = *ts;    

    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    _fluxIn = 0.0;
    _fluxOut = 0.0;

    const T_Scalar electron_effmass = _constants["electron_effmass"];
    const T_Scalar temperature = _constants["OP_temperature"];
    const T_Scalar electron_capture_cross = _constants["electron_capture_cross"];
    //const T_Scalar voltage = _constants["voltage"];
    //const T_Scalar substrate_voltage = _constants["substrate_voltage"];
    //const T_Scalar membrane_voltage = _constants["membrane_voltage"];
    const T_Scalar fermilevelsubstrate = -_constants["substrate_workfunction"] - _constants["substrate_voltage"];
    //const T_Scalar fermilevelmembrane = -_constants["substrate_workfunction"] - _constants["membrane_voltage"];
    //const T_Scalar& dielectric_ionization = _constants["dielectric_ionization"];
    const int subID = _constants["substrate_id"];
    //const int memID = _constants["membrane_id"];
    const int nLevel = _constants["nLevel"];
    const int normal = _constants["normal_direction"];
    const int nTrap = _constants["nTrap"];

    const vector<double> electron_trapdepth = _constants.electron_trapdepth;
    const vector<double> electron_trapdensity = _constants.electron_trapdensity;

    if (int(electron_trapdepth.size()) != nTrap || int(electron_trapdensity.size()) != nTrap)
      throw CException ("trap depth vector size error!");

    T_Scalar fluxCoeff(0), fermilevel(0), scatterfactor(0);
    //T_Scalar sourceTunneling(0);

    for(int c=0; c < cells.getCount(); c++){
      transmission[c] = 0.0;
    }
   
    

    //=======================================//
    // tunneling from substrate to dielectric 
    //=======================================//

    //const T_Scalar energystep = 0.1 * fabs(substrate_voltage - membrane_voltage) / nLevel;
    const T_Scalar energystep = 0.01;
    const T_Scalar alpha = 4.0 * PI * (electron_effmass*ME) / pow(H_SI, 3.0);
    
    fermilevel = fermilevelsubstrate;   

#if DEBUG
    shared_ptr<IntArray> mark(new IntArray(cells.getCount()));
    *mark = 0;
#endif 
   
    for (T_Scalar en = fermilevel-4.0; en <= fermilevel+4.0; en += energystep){

      const T_Scalar supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);

      const T_Scalar fermifunction = FermiFunction(en, fermilevel, temperature);

      //========= transmission coefficient calculation ==========//
      // this scheme only works for Cartesian mesh 

      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	  const FaceGroup& fg = *fgPtr;
	  
	  if (fg.id == subID){ 
	  
	    const StorageSite& faces = fg.site;	 
	    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	    const CRConnectivity& cellCells = mesh.getCellCells();
	    const int nFaces = faces.getCount();
	    
	    for(int f=0; f<nFaces; f++){

	      int c0 = faceCells(f,0);         //the first cell adjacent to boundary
	      int c1 = faceCells(f,1);         //boundary cell

	      transmission[c1] = 1.0;	     
	      
	      int low = c1;
	      int me = c0;	      
	      int high = c0;	 	     
	     
	      for(int l=0; l< nLevel; l++){
		
		T_Scalar dX = cellCentroid[me][normal] - cellCentroid[low][normal];
		T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
		T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
		//T_Scalar valueLow = PositiveValueOf( conduction_band[low] - en);
		//T_Scalar avg = (valueMe + valueLow) / 2.0;
		T_Scalar avg = valueMe;
		T_Scalar exponent = factor * sqrt(avg) * fabs(dX);
		
#if DEBUG
		(*mark)[me] = 1;
#endif		
		transmission[me] = transmission[low] * exp(exponent);
	
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

		if (neighborUp < nCells) {
		  high = neighborUp;
		  low = me;
		  me = high;
		}	
	      }
	    }
	  }
	}
    
      //========= tunneling  calculation ==========//

      for(int c=0; c<nCells; c++){

	for (int i=0; i<nTrap; i++){

	  const T_Scalar stcap = electron_capture_cross * cellVolume[c]; 

	  const T_Scalar endiff = en - (conduction_band[c]-electron_trapdepth[i]);
	
	  // tunneling from substrate to traps 
	
	  if (en-conduction_band[c] < 0){    //this condition determines tunneling only happens close to contact
	    
	    if (endiff < 0)
	      scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	    else scatterfactor = 1.0;
	 
	    fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * 
	      fermifunction * scatterfactor * energystep * QE;
	  
	    _fluxIn += (fluxCoeff * (electron_trapdensity[i] - xCell[c][i])); 
	    rCell[c][i] += (fluxCoeff * (electron_trapdensity[i] - xCell[c][i])); 
	    diag[c](i,i) -= fluxCoeff;
	    
	  }
	  // tunneling from traps to substrate 

	  if(en - conduction_band[c] < 0){
	    if (endiff > 0)
	      scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	    else scatterfactor = 1.0;
	  
	    fluxCoeff = alpha * stcap *  transmission[c] * supplyfunction * 
	      (1-fermifunction) * scatterfactor * energystep * QE;
	  
	    _fluxOut +=  (fluxCoeff * (- xCell[c][i])); 
	    rCell[c][i] += (fluxCoeff * (- xCell[c][i])); 
	    diag[c](i,i) -= fluxCoeff;
	 
	  }
	}
      }
    }
#if 0
    //=======================================//
    // tunneling from membrane to dielectric 
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
	      
	      int c0 = faceCells(f,0);
	      int c1 = faceCells(f,1);
	      
	      transmission[c1] = 1.0;	     
	      
	      int low = c1;
	      int me = c0;	      
	      int high = c0;	      
	      	     
	      for(int l=0; l<nLevel; l++){
		
		T_Scalar dX = cellCentroid[me][normal] - cellCentroid[low][normal];
		T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
		T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
		T_Scalar valueLow = PositiveValueOf( conduction_band[low] - en);
		T_Scalar avg = (valueMe + valueLow) / 2.0;
		T_Scalar exponent = factor * sqrt(avg) * fabs(dX);
		transmission[me] = transmission[low] * exp(exponent);
	
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
		if (neighborUp < nCells) {
		  high = neighborUp;
		  low = me;
		  me = high;
		}	
	      }
	    }
	  }
	}

      for(int c=0; c<nCells; c++){
	
	for(int i=0; i<nTrap; i++){

	  const T_Scalar stcap = electron_capture_cross * cellVolume[c];  

	  const T_Scalar endiff = en - (conduction_band[c]-electron_trapdepth[i]);
	  
	  // tunneling from membrane to traps 
	
	  if (en-conduction_band[c] < 0){

	    if (endiff < 0)
	      scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	    else scatterfactor = 1.0;
	  
	    fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * 
	      fermifunction * scatterfactor * energystep * QE;
	  
	    rCell[c][i] += (fluxCoeff * (electron_trapdensity[i] - xCell[c][i])); 
	    diag[c](i,i) -= fluxCoeff;
	  }
	  
	  // tunneling from traps to membrane
	  if(en - conduction_band[c] < 0 ){
	    if (endiff > 0)
	      scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	    else scatterfactor = 1.0;
	
	    fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * 
	      (1-fermifunction) * scatterfactor * energystep * QE;
	  
	    rCell[c][i] += (fluxCoeff * (-xCell[c][i])); 
	    diag[c](i,i) -= fluxCoeff;
	  }
	}
      }
    }
    
#endif
    //cout << _fluxIn << "  " << _fluxOut << endl;
  }
 
 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _conductionbandField;
  const ElectricModelConstants<T_Scalar>& _constants;
  T_Scalar& _fluxIn;
  T_Scalar& _fluxOut;
};


#endif
