// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _UNDERRELAXER_H_
#define _UNDERRELAXER_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"

  
template<class X, class Diag, class OffDiag>
class Underrelaxer : public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;


  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;

  Underrelaxer(const MeshList& meshes,
               Field& varField,
               const T_Scalar urf) :
    Discretization(meshes),
    _varField(varField),
    _urf(urf)
  {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField&, MultiField&)
  {
    const StorageSite& cells = mesh.getCells();

    const MultiField::ArrayIndex cIndex(&_varField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cIndex,cIndex));
    DiagArray& diag = matrix.getDiag();
    const int nCells = cells.getSelfCount();

    for(int c=0; c<nCells; c++)
    {
        diag[c] /= _urf;
    }
  }
private:
  const Field& _varField;
  const T_Scalar _urf;
};

#endif
