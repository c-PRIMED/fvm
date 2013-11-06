// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLUXJACOBIANMATRIX_H_
#define _FLUXJACOBIANMATRIX_H_

#include "Matrix.h"
#include "CRConnectivity.h"
#include "Array.h"
#include "StorageSite.h"

template<class OffDiag, class X>
class FluxJacobianMatrix : public Matrix
{
public:
  typedef Array<OffDiag> OffDiagArray;
  typedef Array<X> XArray;
  
  FluxJacobianMatrix(const CRConnectivity& conn) :
    Matrix(),
    _conn(conn),
    _coeffL(_conn.getRowDim()),
    _coeffR(_conn.getRowDim())
  {}


  virtual ~FluxJacobianMatrix(){}
  
  DEFINE_TYPENAME("FluxJacobianMatrix<"
             +NumTypeTraits<OffDiag>::getTypeName()+","
             +NumTypeTraits<X>::getTypeName()
             +">");

  virtual void multiply(IContainer& yB, const IContainer& xB) const
  {
    XArray& y = dynamic_cast<XArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    
    const int nRows = _conn.getRowDim();
    for(int nr=0; nr<nRows; nr++)
    {
        const int c0=_conn(nr,0);
        const int c1=_conn(nr,1);
        
        y[nr] = _coeffL[nr]*x[c0] + _coeffR[nr]*x[c1];
    }
  }

  virtual void transpose() {}
  
  virtual void multiplyAndAdd(IContainer& yB, const IContainer& xB) const
  {
    XArray& y = dynamic_cast<XArray&>(yB);
    const XArray& x = dynamic_cast<const XArray&>(xB);
    
    const int nRows = _conn.getRowDim();
    for(int nr=0; nr<nRows; nr++)
    {
        const int c0=_conn(nr,0);
        const int c1=_conn(nr,1);
        
        y[nr] += _coeffL[nr]*x[c0] + _coeffR[nr]*x[c1];
    }
  }

  void setCoeffL(const int f, const OffDiag& c) {_coeffL[f]=c;}
  void setCoeffR(const int f, const OffDiag& c) {_coeffR[f]=c;}
  
  const OffDiag& getCoeffL(const int f) const {return _coeffL[f];}
  const OffDiag& getCoeffR(const int f) const {return _coeffR[f];}
  
  const CRConnectivity& getConnectivity() const {return _conn;}

  virtual void initAssembly()
  {
    _coeffL.zero();
    _coeffR.zero();
  }

private:
  const CRConnectivity& _conn;
  Array<OffDiag> _coeffL;
  Array<OffDiag> _coeffR;
};

#endif
