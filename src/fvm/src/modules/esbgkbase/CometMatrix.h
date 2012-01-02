#ifndef _COMETMATRIX_H_
#define _COMETMATRIX_H_

#include "MatrixJML.h"
#include "Array.h"
#include "ArrayBase.h"
#include "NumType.h"
#include "SquareMatrixESBGK.h"
#include <omp.h>

template<class T>
class CometMatrix : public MatrixJML<T>
{
 public:
  typedef Array<T> TArray;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;

 CometMatrix(const int order):
  _elements(7*order-18),
    _values(_elements),
    _order(order),
    _AMat(3),
    _vel(3)
    {
      _values=0.;
      _AMat.zero();
      _vel=0;
    }

  T& getElement(const int i,const int j)
    {
      int index;
      if(i==j)
      {
	  index=(i-1);
	  return _values[index];
      }
      else if(j==_order)
      {
	  index=_order+3*i-1;
	  return _values[index];
      }
      else if(j==_order-1)
      {
          index=_order+3*i-2;
          return _values[index];
      }
      else if(j==_order-2)
      {
          index=_order+3*i-3;
          return _values[index];
      }
      else if(i==_order)
      {
	  index=6*_order+j-16;
	  return _values[index];
      }
      else if(i==_order-1)
      {
          index=5*_order+j-13;
          return _values[index];
      }
      else if(i==_order-2)
      {
          index=4*_order+j-10;
          return _values[index];
      }
      else
      {
	  throw CException("Invalid index for Comet matrix");
	  return _values[0];
      }
    }

  void printElement(const int& i,const int& j)
    {
      int index;
      if(i==j)
	{
	  index=(i-1);
	  cout<<_values[index]<<endl;
	}
      else if(j==_order)
	{
	  index=i-1+_order;
	  cout<<_values[index]<<endl;
	}
      else if(i==_order)
	{
	  index=2*_order+j-2;
	  cout<<_values[index]<<endl;
	}
      else
	{
	  throw CException("Invalid index for Comet matrix");
	}
    }

  void Solve(TArray& bVec)
  {

    //Replaces bVec with solution vector.
    
    T an1i;
    T an2i;
    T an3i;
    T ain1;
    T ain2;
    T ain3;
    T aii;
    T alpha0;
    T alpha1;
    T alpha2;

    alpha0=0.;
    alpha1=0.;
    alpha2=0.;

    //#pragma omp parallel for default(shared) private(i,an1i,an2i,an3i,aii,ain1,ain2,ain3)
    {
      for(int i=1;i<_order-2;i++)
      {
	  an1i=getElement(_order-2,i);
	  an2i=getElement(_order-1,i);
	  an3i=getElement(_order,i);
	  aii=getElement(i,i);
	  ain1=getElement(i,_order-2);
	  ain2=getElement(i,_order-1);
	  ain3=getElement(i,_order);

	  _AMat.getElement(1,1)-=an1i*ain1/aii;
	  _AMat.getElement(1,2)-=an1i*ain2/aii;
	  _AMat.getElement(1,3)-=an1i*ain3/aii;

	  _AMat.getElement(2,1)-=an2i*ain1/aii;
	  _AMat.getElement(2,2)-=an2i*ain2/aii;
	  _AMat.getElement(2,3)-=an2i*ain3/aii;

	  _AMat.getElement(3,1)-=an3i*ain1/aii;
	  _AMat.getElement(3,2)-=an3i*ain2/aii;
	  _AMat.getElement(3,3)-=an3i*ain3/aii;

	  alpha0+=an1i*bVec[i-1]/aii;
	  alpha1+=an2i*bVec[i-1]/aii;
	  alpha2+=an3i*bVec[i-1]/aii;
      }
    }

    _AMat.getElement(1,1)+=getElement(_order-2,_order-2);
    _AMat.getElement(2,2)+=getElement(_order-1,_order-1);
    _AMat.getElement(3,3)+=getElement(_order,_order);
    
    _vel[0]=(bVec[_order-3]-alpha0);
    _vel[1]=(bVec[_order-2]-alpha1);
    _vel[2]=(bVec[_order-1]-alpha2);

    _AMat.Solve(_vel);
    
    bVec[_order-3]=_vel[0];
    bVec[_order-2]=_vel[1];
    bVec[_order-1]=_vel[2];

    //#pragma omp parallel for default(shared) private(i,aii,ain1,ain2,ain3)
    {
      for(int i=1;i<_order-2;i++)
      {
	  ain1=getElement(i,_order-2);
	  ain2=getElement(i,_order-1);
	  ain3=getElement(i,_order);  
	  aii=getElement(i,i);
	  bVec[i-1]=(bVec[i-1]-ain1*bVec[_order-3]-ain2*bVec[_order-2]-ain3*bVec[_order-1])/aii;
      }
    }
  }

  void zero()
  {
    for(int i=0;i<_elements;i++)
      _values[i]=T_Scalar(0);
  }

  T getTraceAbs()
  {
    T trace=0.;
    for(int i=0;i<_order;i++)
      trace+=fabs(_values[i]);
    return trace;
  }

 private:
  int _elements;
  TArray _values;
  int _order;
  SquareMatrixESBGK<T> _AMat;
  TArray _vel;

};

#endif
