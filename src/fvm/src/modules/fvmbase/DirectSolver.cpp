// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


#include <umfpack.h>

#include "DirectSolver.h"
#include "CRMatrix.h"

DirectSolver::DirectSolver()
{}

DirectSolver::~DirectSolver()
{}

void
DirectSolver::cleanup()
{}


MFRPtr
DirectSolver::solve(LinearSystem & ls)
{

  typedef CRMatrix<double,double,double>  FlatMatrix;
  
  const MultiFieldMatrix& mfMatrix = ls.getMatrix();
  const MultiField& bField = ls.getB();
  MultiField& deltaField = ls.getDelta();

  if (bField.getLength() != 1)
    throw CException("direct solver only works for single matrix");

  const MultiField::ArrayIndex rowIndex = bField.getArrayIndex(0);

  const Matrix& origMatrix = mfMatrix.getMatrix(rowIndex,rowIndex);
  const CRConnectivity& origConn = origMatrix.getConnectivity();

  const ArrayBase& origB = bField[rowIndex];
  const int blockSize = origB.getDataSize()/(sizeof(double)*origB.getLength());

  shared_ptr<CRConnectivity> flatConnPtr = origConn.getMultiTranspose(blockSize);

  FlatMatrix flatMatrix(*flatConnPtr);

  origMatrix.setFlatMatrix(flatMatrix);

  Array<double>& flatCoeffs = flatMatrix.getOffDiag();
  const Array<int>& flatRow = flatConnPtr->getRow();
  const Array<int>& flatCol = flatConnPtr->getCol();
  const ArrayBase& flatB = origB;
  
  ArrayBase& flatDelta = deltaField[rowIndex];
  
  const int nFlatEqs = flatConnPtr->getRowDim();
    
  // original system is in delta form

  mfMatrix.computeResidual(deltaField,bField,ls.getResidual());

  MFRPtr rNorm0(ls.getResidual().getOneNorm());

  if (verbosity >0)
    cout << 0 << ": " << *rNorm0 << endl;
  
  
  double *null = (double *) NULL ;
  
  void *Symbolic, *Numeric ;
  (void) umfpack_di_symbolic (nFlatEqs, nFlatEqs,
                              (int*)flatRow.getData(),
                              (int*)flatCol.getData(),
                              (double*) flatCoeffs.getData(),
                              &Symbolic, null, null) ;
  (void) umfpack_di_numeric ( (int*)flatRow.getData(),
                              (int*)flatCol.getData(),
                              (double*) flatCoeffs.getData(),
                              
                              Symbolic, &Numeric, null, null) ;
  umfpack_di_free_symbolic (&Symbolic) ;
  (void) umfpack_di_solve (UMFPACK_A,
                           (int*)flatRow.getData(),
                           (int*)flatCol.getData(),
                           (double*) flatCoeffs.getData(),
                           (double*) flatDelta.getData(),
                           (double*) flatB.getData(),
                           Numeric, null, null) ;
   umfpack_di_free_numeric (&Numeric) ;

   // umfpack solves ax=b, we want ax+b=0;
   double *flatDeltaData = (double*) flatDelta.getData();
   for(int n=0; n<nFlatEqs; n++)
     flatDeltaData[n] *= -1.0;
   
  mfMatrix.computeResidual(deltaField,bField,ls.getResidual());

  MFRPtr rNormN(ls.getResidual().getOneNorm());

  if (verbosity >0)
    cout <<  "Final : " << *rNormN << endl;
  

  return rNorm0;
}

void
DirectSolver::smooth(LinearSystem& ls)
{
  throw CException("cannot use DirectSolver as preconditioner");
}
