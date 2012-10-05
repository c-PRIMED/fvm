// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ONEDCONDUCTION_H_
#define _ONEDCONDUCTION_H_

#include "Array.h"

template<class T>
void TDMA(Array<T>& ap, Array<T>& ae, Array<T>& aw, Array<T>& b, Array<T>& x)
{
  const int n=ap.getLength();
  for(int i=1; i<n; i++)
  {
      const T m(aw[i]/ap[i-1]);
      ap[i] -= m*ae[i-1];
      b[i] -= m*b[i-1];
  }
  x[n-1] = -b[n-1]/ap[n-1];
  for(int i=n-2; i>=0; i--)
    x[i] = (b[i] - ae[i]*x[i+1])/ap[i];
}

template<class T>
class OneDConduction
{
public:
  OneDConduction(const int nCells, const T& kConst) :
    _nCells(nCells),
    _kConst(kConst),
    _xL(0),
    _xR(1),
    _x(new Array<T>(_nCells))
  {}
  
  void solve()
  {
    Array<T> ae(_nCells);
    Array<T> aw(_nCells);
    Array<T> ap(_nCells);
    Array<T> b(_nCells);

    ae = T(0.);
    aw = T(0.);
    ap = T(0.);
    b = T(0.);
    
    const T dx = T(1.0)/_nCells;
    
    for(int nf=0; nf<=_nCells; nf++)
    {
        int c0 = nf-1;
        int c1 = nf;

        const T xf =  dx*nf;
        const T kf = 1.0 + _kConst*xf;
        T coeff = kf/dx;
        if (nf==0)
        {
            coeff *= 2.;
            c0 = 0;
            b[c0] += coeff*_xL;
            ap[c0] -= coeff;
        }
        else if (nf==_nCells)
        {
            coeff *=2;
            b[c0] += coeff*_xR;
            ap[c0] -= coeff;
        }
        else
        {
            ae[c0] = coeff;
            aw[c1] = coeff;
            ap[c0] -= coeff;
            ap[c1] -= coeff;
        }
    }

    cout << " ae[0] " << ae[0] << endl;
    TDMA(ap,ae,aw,b,*_x);
  }

  boost::shared_ptr<Array<T> > getSolution() {return _x;}
  
private:
  const int _nCells;
  T _kConst;
  T _xL;
  T _xR;
  boost::shared_ptr<Array<T> > _x;
  
};
#endif
