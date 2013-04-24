// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BATTERYPC_BCS_H_
#define _BATTERYPC_BCS_H_

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
class BatteryPC_BCS
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

 
 BatteryPC_BCS(const StorageSite& faces,
             const Mesh& mesh,
             const GeomFields& geomFields,
             Field& varField,
             Field& fluxField,
             MultiFieldMatrix& matrix,
             MultiField& xField, MultiField& rField) :
    _faces(faces),
    _cells(mesh.getCells()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])),
    _faceCells(mesh.getFaceCells(_faces)),
    _varField(varField),
    _fluxField(fluxField),
    _xIndex(&_varField,&_cells),
    _fluxIndex(&_fluxField,&_faces),
    _dRdX(dynamic_cast<CCMatrix&>(matrix.getMatrix(_xIndex,_xIndex))),
    _dFluxdX(dynamic_cast<FMatrix&>(matrix.getMatrix(_fluxIndex,_xIndex))),
    _dFluxdFlux(dynamic_cast<BBMatrix&>(matrix.getMatrix(_fluxIndex,_fluxIndex))),
    _assembler(_dRdX.getPairWiseAssembler(_faceCells)),
    _dRdXDiag(_dRdX.getDiag()),
    _x(dynamic_cast<XArray&>(xField[_xIndex])),
    _r(dynamic_cast<XArray&>(rField[_xIndex])),
    _flux(dynamic_cast<XArray&>(xField[_fluxIndex])),
    _rFlux(dynamic_cast<XArray&>(rField[_fluxIndex])),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _is2D(mesh.getDimension()==2)    
  {}
    
    void applySingleEquationDirichletBC(int f, const T_Scalar& bValue, const int v) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    // the current value of flux and its Jacobians
    const T_Scalar fluxB = -(_r[c1])[v];
    //const T_Scalar dFluxdXC0 = -(_assembler.getCoeff10(f))(v,v);
    const OffDiag dFluxdXC0_orig = -_assembler.getCoeff10(f);
    const T_Scalar dFluxdXC1 = -(_dRdXDiag[c1])[v];
    OffDiag dRC0dXC1_orig = _assembler.getCoeff01(f);
    const T_Scalar dRC0dXC1 = dRC0dXC1_orig[v];
    
    // since we know the boundary value, compute the boundary
    // x correction and it's contribution to the residual for c0; we
    // can then eliminate the coefficient to the boundary cell
    
    const T_Scalar dXC1 = bValue - (_x[c1])[v];
    const T_Scalar dFlux = dFluxdXC1*dXC1;
    const T_Scalar dRC0 = dRC0dXC1*dXC1;
    (_r[c0])[v] += dRC0;  
    
    (_assembler.getCoeff01(f))[v] = NumTypeTraits<T_Scalar>::getZero();

    // set the boundary value and make its correction equation an
    // identity
    (_x[c1])[v] = bValue;
    (_assembler.getCoeff10(f))[v] = NumTypeTraits<T_Scalar>::getZero();
    (_r[c1])[v] = NumTypeTraits<T_Scalar>::getZero();
    (_dRdXDiag[c1])[v] = NumTypeTraits<T_Scalar>::getNegativeUnity();

    //setup the equation for the boundary flux correction
    _dFluxdX.setCoeffL(f,dFluxdXC0_orig);
    dRC0dXC1_orig[v] = NumTypeTraits<T_Scalar>::getZero();
    _dFluxdX.setCoeffR(f,dRC0dXC1_orig);
    (_flux[f])[v] = fluxB;
    (_rFlux[f])[v] = dFlux;
    Diag _dFluxdFlux_orig = _dFluxdFlux[f];
    _dFluxdFlux_orig[v] = NumTypeTraits<T_Scalar>::getNegativeUnity();
    _dFluxdFlux[f] = _dFluxdFlux_orig;
  }

  void applySingleEquationDirichletBC(const FloatValEvaluator<T_Scalar>& bValue, const int v) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applySingleEquationDirichletBC(i,bValue[i],v);
  }
  
  void applySingleEquationNeumannBC(const int f,
                      const T_Scalar& specifiedFlux, const int v) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    // the current value of flux and its Jacobians
    const T_Scalar fluxB = -(_r[c1])[v];
        
        
    // since we know the boundary flux, compute the boundary flux
    // correction and add it to the c0 residual; also eliminate
    // coeff to the boundary cell and remove it from the ap coeff
    
    const T_Scalar dFlux = specifiedFlux*_faceAreaMag[f] - fluxB;
    // setup the equation for the boundary value; the coefficients
    // are already computed so just need to set the rhs
    (_r[c1])[v] = dFlux;

    // mark this row as a "boundary" row so that we will update it
    // after the overall system is solved
    _dRdX.setBoundary(c1);

    // fix the boundary flux to the specified value
    (_flux[f])[v] = specifiedFlux*_faceAreaMag[f];

    Diag _dFluxdFlux_orig = _dFluxdFlux[f];
    _dFluxdFlux_orig[v] = NumTypeTraits<T_Scalar>::getNegativeUnity();
    _dFluxdFlux[f] = _dFluxdFlux_orig;
  }
  
  void applySingleEquationNeumannBC(const T_Scalar& bFlux, const int v) const
  {
    for(int i=0; i<_faces.getCount(); i++)
      applySingleEquationNeumannBC(i,bFlux,v);
  }
  
protected:
  const StorageSite& _faces;
  const StorageSite& _cells;
  const IntArray& _ibType;
  const CRConnectivity& _faceCells;
  const Field& _varField;
  const Field& _fluxField;
  const MultiField::ArrayIndex _xIndex;
  const MultiField::ArrayIndex _fluxIndex;
  CCMatrix& _dRdX;
  FMatrix& _dFluxdX;
  BBMatrix& _dFluxdFlux;
  CCAssembler& _assembler;
  DiagArray& _dRdXDiag;
  XArray& _x;
  XArray& _r;
  XArray& _flux;
  XArray& _rFlux;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  const bool _is2D;
  
};

#endif

