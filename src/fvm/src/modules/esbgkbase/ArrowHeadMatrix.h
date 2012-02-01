#ifndef _ARROWHEADMATRIX_H_
#define _ARROWHEADMATRIX_H_

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "NumType.h"
#include "Array.h"
#include "ArrayBase.h"
#include "Vector.h"
#include "VectorTranspose.h"
#include "SquareTensor.h"
#include "StorageSite.h"

template<class X, int K>
class ArrowHeadMatrix
{
public:
  
  typedef Array<X> XArray;
  typedef Array<Vector<X,K> > VectorXKArray;
  typedef Array<VectorTranspose<X,K> > TVectorXKArray;
  typedef SquareTensor<X,K> TensorXK;
  typedef Vector<X,K> VectorXK;
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  ArrowHeadMatrix(const int order) :
    _order(order),
    _numDir(order-3),
    _d(_numDir),
    _r(_numDir),
    _c(_numDir),
    _l(3),
    _bl(),
    _xl()
  {
    _d.zero();
    _r.zero();
    _c.zero();
    _l.zero();
    _bl.zero();
    _xl.zero();
  }

  T_Scalar& getElement(const int i,const int j)
  {
      int index;	
      if(i<=_numDir)
      {
	  index=(i-1);
	  if(i==j)
	      return _d[index];
	  else
	      return _c[index][j-_order+2];
      }
      else
      { 
	  index=(j-1);
	  if(j<=_order-3)
	      return (_r[index])[i-_order+2];
	  else
	    return _l(i-_order+2,j-_order+2);
      }
  }

  void Solve(XArray& bVec)
  {
    _bl[0]=bVec[_order-3];
    _bl[1]=bVec[_order-2];
    _bl[2]=bVec[_order-1];
    for(int i=0; i<_numDir; i++)
    {
	_l-=((_c[i]).getTensor(_r[i]))/_d[i];
	_bl-=_r[i]*bVec[i]/_d[i];
    }

    _xl=_bl/_l;
 
    for(int i=0; i<_order-3; i++)
      bVec[i]=(bVec[i]-_c[i]*_xl)/_d[i];
    bVec[_order-3]=_xl[0];
    bVec[_order-2]=_xl[1];
    bVec[_order-1]=_xl[2];        
  }

  void zero()
  {
    _d.zero();
    _r.zero();
    _c.zero();
    _l.zero();
  }

  T_Scalar getTraceAbs()
  {
    T_Scalar trace=0.;
    for(int i=0;i<_order-3;i++)
      trace+=fabs(_d[i]);
    trace+=(fabs(_l(0,0))+fabs(_l(1,1))+fabs(_l(2,2)));
    return trace;
  }

  
private:

  const int _order;
  const int _numDir;
  XArray _d;
  VectorXKArray _r;
  TVectorXKArray _c;
  TensorXK _l;
  VectorXK _bl;
  VectorXK _xl;
};

#endif
