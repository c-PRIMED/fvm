// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GENERICPHONONBCS_H_
#define _GENERICPHONONBCS_H_

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
class BaseGenericPhononBCS
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<T_Scalar> TArray;
  typedef Array<int> IntArray;
  
  typedef Vector<T_Scalar,3> VectorT3;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef FluxJacobianMatrix<Diag,X> FMatrix;
  typedef DiagonalMatrix<Diag,X> BBMatrix;

  typedef Array<Diag> DiagArray;
  typedef Array<OffDiag> OffDiagArray;
  
  typedef Array<X> XArray;
  typedef Array<VectorT3> VectorT3Array;
  
  
 BaseGenericPhononBCS(const StorageSite& faces,
		      const Mesh& mesh,
		      const GeomFields& geomFields,
		      Field& varField,
		      MultiFieldMatrix& matrix,
		      MultiField& xField, MultiField& rField) :
  _faces(faces),
    _cells(mesh.getCells()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])),
    _faceCells(mesh.getFaceCells(_faces)),
    _varField(varField),
    _xIndex(&_varField,&_cells),
    _dRdX(dynamic_cast<CCMatrix&>(matrix.getMatrix(_xIndex,_xIndex))),
    _assembler(_dRdX.getPairWiseAssembler(_faceCells)),
    _dRdXDiag(_dRdX.getDiag()),
    _x(dynamic_cast<XArray&>(xField[_xIndex])),
    _r(dynamic_cast<XArray&>(rField[_xIndex])),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces]))
      {}
  
  void applyDirichletBC(int f, const X& bValue) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const OffDiag dRC0dXC1 = _assembler.getCoeff01(f);
    
    // since we know the boundary value, compute the boundary
    // x correction and it's contribution to the residual for c0; we
    // can then eliminate the coefficient to the boundary cell
    
    const X dXC1 = bValue - _x[c1];
    const X dRC0 = dRC0dXC1*dXC1;
    _r[c0] += dRC0;
    
    _assembler.getCoeff01(f) = NumTypeTraits<OffDiag>::getZero();

    // set the boundary value and make its correction equation an
    // identity
    _x[c1] = bValue;
    _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getZero();
    _r[c1] = NumTypeTraits<X>::getZero();
    _dRdXDiag[c1] = NumTypeTraits<Diag>::getNegativeUnity();

 }

  void applyDirichletBC(const X& bValue) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDirichletBC(i,bValue);
  }
  
  void applyDirichletBC(const FloatValEvaluator<X>& bValue) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyDirichletBC(i,bValue[i]);
  }
  

  void applyExtrapolationBC() const
  {
    for(int i=0; i<_faces.getCount(); i++)
	applyExtrapolationBC(i);
  }

  // boundary value = cell value, flux as defined by interior discretization
  
  void applyExtrapolationBC(const int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
        
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    // the current value of flux and its Jacobians
    const Diag dFluxdXC1 = -_dRdXDiag[c1];

    const X xc0mxc1 = _x[c0]-_x[c1];
        
    // eliminate boundary dependency from cell equation
    _dRdXDiag[c0] += dFluxdXC1;
    _r[c0] += dFluxdXC1*xc0mxc1;
    _assembler.getCoeff01(f) = 0;

    // boundary value equation
    _dRdXDiag[c1] = NumTypeTraits<Diag>::getNegativeUnity();
    _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getUnity();
    _r[c1] = xc0mxc1;
    _dRdX.setBoundary(c1);
  
 }
  

  void applyInterfaceBC() const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyInterfaceBC(i);
  }

  void applyFlowBC(const TArray& convFlux, const X& bValue) const
  {
    for(int f=0; f<_faces.getCount(); f++)
      if (convFlux[f] < 0)
        applyDirichletBC(f,bValue);
      else
        applyExtrapolationBC(f);
  }

  void applyNonzeroDiagBC() const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applyNonzeroDiagBC(i);
  }

  void applyNonzeroDiagBC(int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;    
   
    _dRdXDiag[c1][0] = T_Scalar(-1.0);   
  }

  
protected:
  const StorageSite& _faces;
  const StorageSite& _cells;
  const IntArray& _ibType;
  const CRConnectivity& _faceCells;
  const Field& _varField;
  const MultiField::ArrayIndex _xIndex;
  CCMatrix& _dRdX;
  CCAssembler& _assembler;
  DiagArray& _dRdXDiag;
  XArray& _x;
  XArray& _r;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
};

#endif

