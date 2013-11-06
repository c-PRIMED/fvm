// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _DIAGONALMATRIX_H_
#define _DIAGONALMATRIX_H_



#include "Matrix.h"
#include "CRConnectivity.h"
#include "Array.h"


template<class Diag, class X>
class DiagonalMatrix : public Matrix
{
public:
  typedef Array<Diag> DiagArray;
  typedef Array<X> XArray;
  
  DiagonalMatrix(const int length) :
    Matrix(),
    _length(length),
    _diag(_length)
  {
    logCtor();
  }

  virtual ~DiagonalMatrix(){}
  
  DEFINE_TYPENAME("DiagonalMatrix<"
                  +NumTypeTraits<Diag>::getTypeName()+","
                  +NumTypeTraits<X>::getTypeName()
                  +">");
  
  virtual void multiply(IContainer& yB, const IContainer& xB) const
  {
    XArray& y = dynamic_cast<XArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    for(int i=0; i<_length; i++)
    {
        y[i] = _diag[i]*x[i];
    }
  }

  virtual void multiplyAndAdd(IContainer& yB, const IContainer& xB) const
  {
    XArray& y = dynamic_cast<XArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    for(int i=0; i<_length; i++)
    {
        y[i] += _diag[i]*x[i];
    }
  }
  
  virtual void forwardGS(IContainer& xB, IContainer& bB, IContainer&) const
  {
    XArray& x = dynamic_cast<XArray&>(xB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    for(int i=0; i<_length; i++)
      x[i] = -b[i]/_diag[i];
  }
  
  virtual void reverseGS(IContainer& xB, IContainer& bB,IContainer& r) const
  {
    forwardGS(xB,bB,r);
  }

  virtual void solveBoundary(IContainer& xB, IContainer& bB,IContainer& r) const
  {
    forwardGS(xB,bB,r);
  }

  virtual void computeResidual(const IContainer& xB, const IContainer& bB,
                               IContainer& rB) const
  {
    const XArray& x = dynamic_cast<const XArray&>(xB);
    const XArray& b = dynamic_cast<const XArray&>(bB);
    XArray& r = dynamic_cast<XArray&>(rB);
    
    for(int nr=0; nr<_length; nr++)
    {
        r[nr] = b[nr] + _diag[nr]*x[nr];
    }
  }

  virtual void transpose() {}
  
  Diag& operator[](int n) {return _diag[n];}
  const Diag& operator[](int n) const {return _diag[n];}

  void addToDiag(const int i, const Diag& d)
  {
    _diag[i] += d;
  }

  void unitize()
  {
    _diag = NumTypeTraits<Diag>::getNegativeUnity();
  }

  void unitize(const int i)
  {
    _diag[i] = NumTypeTraits<Diag>::getNegativeUnity();
  }

  virtual shared_ptr<CRConnectivity>
  createCoarseConnectivity(const IContainer& gCoarseIndex,
                           const CRConnectivity& coarseToFine,
                           const StorageSite& coarseRowSite,
                           const StorageSite& coarseColSite)
  {
    // we don't really need one but all matrices are supposed to supply one 
    return shared_ptr<CRConnectivity>(new CRConnectivity(coarseRowSite,coarseColSite));
  }

  virtual shared_ptr<Matrix>
  createCoarseMatrix(const IContainer& gCoarseIndex,
                     const CRConnectivity& coarseToFine,
                     const CRConnectivity& )
  {
    const int nCoarseRows = coarseToFine.getRowDim();

    shared_ptr<DiagonalMatrix> coarseMatrix(new DiagonalMatrix(nCoarseRows));

    Array<Diag>& coarseDiag = coarseMatrix->_diag;
    coarseDiag.zero();

    for(int nrCoarse=0; nrCoarse<nCoarseRows; nrCoarse++)
    {
        // loop over the fine rows that make up this coarse row
        const int nFine = coarseToFine.getCount(nrCoarse);
        for(int nfr=0; nfr<nFine; nfr++)
        {
            const int nrFine = coarseToFine(nrCoarse,nfr);

            coarseDiag[nrCoarse] += _diag[nrFine];
            
        }
    }

    return coarseMatrix;
  }

  virtual void initAssembly()
  {
    _diag.zero();
  }

  virtual bool isInvertible() {return true;}

private:
  const int _length;
  Array<Diag> _diag;
};
#endif

