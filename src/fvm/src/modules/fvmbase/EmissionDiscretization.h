// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _EMISSIONDISCRETIZATION_H_
#define _EMISSIONDISCRETIZATION_H_

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

Emission model modifies d00,  x0, x1

*************************************************/


template <class X, class Diag, class OffDiag>
class EmissionDiscretization : public Discretization
{
 public:
  
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef Array<Vector<T_Scalar,3> > VectorT3Array;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;
  
  EmissionDiscretization(const MeshList& meshes,
			 const GeomFields& geomFields,
			 const Field& varField,
			 const Field& electricField,
			 const ElectricModelConstants<T_Scalar>& constants):
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _electricField(electricField),
    _constants(constants)
    {}


  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const int nCells = cells.getSelfCount();
    
    const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array&> (_electricField[cells]);
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(_varField[cells]);
    
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    //OffDiagArray& offdiag = matrix.getOffDiag();
    
    const T_Scalar optical_dielectric_constant = _constants["optical_dielectric_constant"];
    const T_Scalar temperature = _constants["OP_temperature"];
    const T_Scalar poole_frenkel_emission_frequency = _constants["poole_frenkel_emission_frequency"];
    const int nTrap = _constants["nTrap"];
    vector<T_Scalar> electron_trapdepth = _constants.electron_trapdepth;
    const T_Scalar beta = sqrt( QE / (PI * E0_SI * optical_dielectric_constant) );

    for(int c=0; c<nCells; c++){
      
      for (int i=0; i<nTrap; i++){
	
	T_Scalar expt = (electron_trapdepth[i] - beta * sqrt(mag(electric_field[c]))) * QE / (K_SI * temperature);

	//if (expt < 0.0)
	//throw CException("exponential error in Emission model");
	if (expt > 0.0)
	{

	    T_Scalar fluxCoeff = cellVolume[c] * poole_frenkel_emission_frequency * exp(-expt);
	
	    rCell[c][i] -= (fluxCoeff * xCell[c][i]);

	    diag[c](i,i) -= fluxCoeff;
	    //diag[c][i] -= fluxCoeff;

	    rCell[c][nTrap] += fluxCoeff * xCell[c][i];
	}
      }
    }
  }

 private:
    const GeomFields& _geomFields;
    const Field& _varField;
    const Field& _electricField;
    const ElectricModelConstants<T_Scalar>& _constants;

};


#endif
