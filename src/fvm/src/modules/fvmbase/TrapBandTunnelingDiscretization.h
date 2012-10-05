// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _TRAPBANDTUNNELINGDISCRETIZATION_H_
#define _TRAPBANDTUNNELINGDISCRETIZATION_H_

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

template <class X, class Diag, class OffDiag>
class TrapBandTunnelingDiscretization : public Discretization
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;
  typedef Array<Vector<T_Scalar, 3> > VectorT3Array;

  TrapBandTunnelingDiscretization(const MeshList& meshes,
				  const GeomFields& geomFields,
				  const Field& varField,
				  const Field& electricField,
				  const Field& conductionbandField,
				  const ElectricModelConstants<T_Scalar>& constants):
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
    
    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array&> (_electricField[cells]);

    const VectorT3Array& cellCentroid = 
      dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(_varField[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const CRConnectivity& cellCells = mesh.getCellCells();
    
    TArray* ts = (new TArray(cells.getCount()));
    *ts = 0;
    TArray& transmission = *ts;    
    
    const T_Scalar electron_capture_cross = _constants["electron_capture_cross"];
    const T_Scalar electron_effmass = _constants["electron_effmass"];
    const int normal = _constants["normal_direction"];
    const int nTrap = _constants["nTrap"];
    vector<T_Scalar> electron_trapdepth = _constants.electron_trapdepth;
    //const T_Scalar epsilon = 1e-18;

    const int nMax = 200;
    
    T_Scalar transmissionHigh(0), transmissionLow(0);

    int idHigh, idLow, high, low, me, count;

    bool foundHigh, foundLow, flag;

    idHigh = idLow = 0;

    foundHigh = foundLow = false;

    for(int c=0; c<cells.getCount(); c++){
      transmission[c] = 0.0;
    }
    
    for(int c=0; c<nCells; c++){   

      for(int i=0; i<nTrap; i++){
	
	T_Scalar en = conduction_band[c] - electron_trapdepth[i];

	transmission[c] = 1.0;
	
	//-------------------------------------------------------------------------//

	high = low = me = c;

	flag = false;        //hit the boundary or not?

	count = 0;

	while(flag == false && count < nMax ){

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
	  else flag = true;
	 
	  T_Scalar dX = cellCentroid[me][normal] - cellCentroid[low][normal];
	  T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	  T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
	  T_Scalar valueLow = PositiveValueOf( conduction_band[low] - en);
	  T_Scalar avg = (valueMe + valueLow) / 2.0;
	  T_Scalar exponent = factor * sqrt(avg) * fabs(dX);
	  
	  transmission[me] = transmission[low] * exp(exponent);

	  if (en - conduction_band[me] >0 ){

	    foundHigh = true;
	    idHigh = me;
	    transmissionHigh = transmission[me];
	    //  cout << "found high " << c <<"  " << me << "  " << transmissionHigh << endl;
	    break;

	  }
	  count ++;
	}	

#if 0 
	//-------------------------------------------------------------------------//

	high = low = me = c;
	
	flag = false;
	
	count = 0;
      
	while(flag == false && count < nMax){
	
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
	  else flag = true;

	  T_Scalar dX = cellCentroid[me][normal] - cellCentroid[low][normal];
	  T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	  T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
	  T_Scalar valueLow = PositiveValueOf( conduction_band[low] - en);
	  T_Scalar avg = (valueMe + valueLow) / 2.0;
	  T_Scalar exponent = factor * sqrt(avg) * fabs(dX);
	  
	  transmission[me] = transmission[low] * exp(exponent);
	
	  if (en - conduction_band[me] >0 ){
	    
	    foundLow = true;
	    idLow = me;
	    transmissionLow = transmission[me];
	    //cout << "found low " << c <<"  " << me << endl;
	    break;
	  }
	  //	cout << c << " low " << low << endl; 
	  count ++;
	  
	}
#endif
	//-------------------------------------------------------------------------//

	const T_Scalar ef = mag(electric_field[c]);
	
	const T_Scalar alpha = cellVolume[c] * QE * ef * ef * electron_capture_cross /    
	  (16 * PI*PI * HBAR_SI * electron_effmass * electron_trapdepth[i]);
      
	if (foundLow == true){	
	  rCell[c][i] -= alpha * transmissionLow  * xCell[c][i];
	  diag[c](i,i) -= alpha  * transmissionLow;
	  //diag[c][i] -= alpha  * transmissionLow;
	  rCell[idLow][nTrap] += alpha * transmissionLow * xCell[c][i];
	}

	if (foundHigh == true){
	  rCell[c][i] -= alpha * transmissionHigh * xCell[c][i];
	  diag[c](i,i) -= alpha  *  transmissionHigh;
	  //diag[c][i] -= alpha  *  transmissionHigh;
	  rCell[idHigh][nTrap] += alpha * transmissionHigh * xCell[c][i];
	}
      }
    }
  }	    
	

 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _electricField;
  const Field& _conductionbandField;
  const ElectricModelConstants<T_Scalar>& _constants;
 
	      
};

#endif
