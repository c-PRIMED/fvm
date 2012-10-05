// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LINEARIZEDIELECTRIC_H_
#define _LINEARIZEDIELECTRIC_H_ 

#include "Mesh.h"
#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"




template<class X, class Diag, class OffDiag>
class LinearizeDielectric
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef Array<T_Scalar> TArray;
  typedef Array<X> XArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
 

 LinearizeDielectric( const GeomFields& geomFields,
		      const Field& dielectric_constant,
		      const T_Scalar dielectric_thickness,
		      Field& varField,
		      const T_Scalar source=0.0):
    _geomFields(geomFields),
    _varField(varField), 
    _dielectric_constant(dielectric_constant),
    _dielectric_thickness(dielectric_thickness),
    _source(source)
    {}


  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getParentFaceGroupSite();
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);

    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex)); 
    const TArray& faceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
    
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);

    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const CRConnectivity& cellCells = mesh.getCellCells();

    const TArray& dcCell =
      dynamic_cast<const TArray&>(_dielectric_constant[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    DiagArray& diag = matrix.getDiag();

    const int nCells = cells.getSelfCount();

    for (int c=0; c<nCells; c++)
      {
	const int c0 = c;
	const int nb = cellCells.getCount(c);
	if (nb!=2) 
	  throw CException("invalid connectivity in shellMesh");
	
	//source term
	const T_Scalar volume = faceAreaMag[c0] * _dielectric_thickness;

	const T_Scalar src = 0.0;
	
	rCell[c0] += src*volume;

	//diffusion term
	for (int i=0; i<nb; i++){
	
	  const int c1 = cellCells(c,i);

	  OffDiag& offdiag = matrix.getCoeff(c0,  c1);

	  VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
	  T_Scalar dsMag = mag(ds);
	  
	  const T_Scalar dc = harmonicAverage(dcCell[c0], dcCell[c1]);
	
	  const T_Scalar diffCoeff = dc*faceAreaMag[c0]/(dsMag + 0.5*_dielectric_thickness);

	  const X dFlux = diffCoeff*(xCell[c1]-xCell[c0]);

	  rCell[c0] += dFlux;

	  offdiag = diffCoeff;

	  diag[c0] -= diffCoeff;

	} 
      }
  }
private:
  
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _dielectric_constant;
  const T_Scalar _dielectric_thickness;
  const T_Scalar _source;
  
};

#endif
