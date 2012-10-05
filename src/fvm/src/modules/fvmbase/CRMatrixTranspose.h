// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CRMATRIXTRANSPOSE_H_
#define _CRMATRIXTRANSPOSE_H_


#include "Matrix.h"
#include "CRConnectivity.h"
#include "Array.h"
#include "StorageSite.h"


/**
 * 
 * This class is similar to CRMatrix except that it stores the
 * transpose of the matrix implied by the CRConnectivity object that
 * is given to it. This is useful in situations like node to cell
 * interpolation matrix where the connectivity we have is cell to
 * nodes but the matrix we want is one which stores the weights of
 * node to cell interpolation. In such cases we can create a matrix of
 * this class using the cell to node connectivity but use it in a
 * multiplication operation with the cell value array to get the node
 * value array. It is also very useful for ib faces since storing ib
 * face to cell (or ib face to particle) connectivity is a lot cheaper
 * than the transpose.
 * 
 *
 * The class is templated with the types of coefficient, x and b
 * arrays as template parameters. Since the primary use of this class
 * is for storing interpolation matrices from one site to another we
 * have a special constructor that can create a CRMatrixTranspose
 * object of the same Coeff type but different X and B (e.g. to
 * interpolate vector variables). For this reason the coeff is stored
 * as a shared_ptr.
 * 
 */

template<class T_Coeff, class X, class B>
class CRMatrixTranspose : public Matrix
{
public:
  typedef T_Coeff Coeff;
  typedef Array<Coeff> CoeffArray;
  typedef Array<X> XArray;
  typedef Array<B> BArray;
  

  CRMatrixTranspose(const CRConnectivity& conn) :
    Matrix(),
    _conn(conn),
    _row(_conn.getRow()),
    _col(_conn.getCol()),
    _coeffPtr(new CoeffArray(_col.getLength())),
    _coeff(*_coeffPtr)
  {
    logCtor();
  }

  // meant to be used to create CRMatrixTranspose<T,Vec3,Vec3> from
  // CRMatrix<T,T,T>
  
  template<class X2, class B2>
  CRMatrixTranspose(const CRMatrixTranspose<Coeff,X2,B2>& m) :
    Matrix(),
    _conn(m.getConnectivity()),
    _row(_conn.getRow()),
    _col(_conn.getCol()),
    _coeffPtr(m.getCoeffPtr()),
    _coeff(*_coeffPtr)
  {
    logCtor();
  }

  
  DEFINE_TYPENAME("CRMatrixTranspose<"
                  +NumTypeTraits<Coeff>::getTypeName()+","
                  +NumTypeTraits<X>::getTypeName()+","
                  +NumTypeTraits<B>::getTypeName()
             +">");

  virtual void initAssembly()
  {
    _coeff.zero();
  }

  /**
   * y = this * x
   * 
   */

  virtual void multiply(IContainer& yB, const IContainer& xB) const
  {
    yB.zero();
    this->multiplyAndAdd(yB,xB);
  }

  /**
   * y += this * x
   * 
   */

  virtual void multiplyAndAdd(IContainer& yB, const IContainer& xB) const
  {
    BArray& y = dynamic_cast<BArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    
    const int nRows = _conn.getRowSite().getCount();
    for(int nr=0; nr<nRows; nr++)
    {
        for (int nb = _row[nr]; nb<_row[nr+1]; nb++)
        {
            const int j = _col[nb];
            y[nr] += _coeff[nb]*x[j];
        }
    }
  }

  
  const CRConnectivity& getConnectivity() const {return _conn;}

  Array<Coeff>& getCoeff() {return _coeff;}
  const Array<Coeff>& getCoeff() const {return _coeff;}
  shared_ptr<Array<Coeff> > getCoeffPtr() const {return _coeffPtr;}
  
private:
  const CRConnectivity& _conn;
  const Array<int>& _row;
  const Array<int>& _col;
  shared_ptr<Array<Coeff> > _coeffPtr;
  Array<Coeff>& _coeff;
};

#endif

