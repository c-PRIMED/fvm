// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ARROWHEADMATRIX_H_
#define _ARROWHEADMATRIX_H_

#include "MatrixJML.h"
#include "Array.h"
#include "ArrayBase.h"
#include "NumType.h"
#include <math.h>

template<class T>
class ArrowHeadMatrix : public MatrixJML<T>
{
 public:
  typedef Array<T> TArray;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;

 ArrowHeadMatrix(const int order):
  _elements(3*order-2),
    _values(_elements),
    _order(order)
    {
      _values=0.;
    }

  T& operator()(const int i, const int j)
    {
      int index;
      if(i==j)
	{
	  index=(i-1);
	  return _values[index];
	}
      else if(j==_order)
	{
	  index=i-1+_order;
	  return _values[index];
	}
      else if(i==_order)
	{
	  index=2*_order+j-2;
	  return _values[index];
	}
      else
	{
	  throw CException("Invalid index: Arrowhead matrix");
	  return _values[0];
	}
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
	  index=i-1+_order;
	  return _values[index];
	}
      else if(i==_order)
	{
	  index=2*_order+j-2;
	  return _values[index];
	}
      else
	{
	  throw CException("Invalid index: Arrowhead matrix");
	  return _values[0];
	}
    }

  void printElement(const int& i,const int& j)
    {
      int index;
      if(i==j)
	{
	  index=(i-1);
	  cout<<_values[index];
	}
      else if(j==_order)
	{
	  index=i-1+_order;
	  cout<<_values[index];
	}
      else if(i==_order)
	{
	  index=2*_order+j-2;
	  cout<<_values[index];
	}
      else
	{
	  throw CException("Invalid index for Arrowhead matrix");
	}
    }

  void print()
  {
    for(int i=1;i<_order+1;i++)
      {
	for(int j=1;j<_order+1;j++)
	  {
	    if((i==j)||(i==_order)||(j==_order))
	      {
		printElement(i,j);
		cout<<" ";
	      }
	    else
	      cout<<0<<" ";
	  }
	cout<<endl;
      }
  }

  void Solve(TArray& bVec)
  {

    //Replaces bVec with solution vector.
    
    T ani;
    T ain;
    T aii;
    T alpha;
    T beta;

    alpha=0.;
    beta=0.;

 
    for(int i=1;i<_order;i++)
      {
	ani=(*this)(_order,i);
	aii=(*this)(i,i);
	ain=(*this)(i,_order);
	
	alpha+=ani*bVec[i-1]/aii;
	beta+=ani*ain/aii;
      }
    
    T bn;
    bn=(bVec[_order-1]-alpha)/((*this)(_order,_order)-beta);
    
    for(int i=1;i<_order;i++)
      {
	ain=(*this)(i,_order);
	aii=(*this)(i,i);
	bVec[i-1]=(bVec[i-1]-ain*bn)/aii;
      }
    
    bVec[_order-1]=bn;
  }

  void SolveDiag(TArray& bVec)
  {
    for(int i=1;i<_order;i++)
      bVec[i-1]/=(*this)(i,i);
    
    bVec[_order-1]=0;
  }

  void smoothGS(TArray& bVec)
  {
    T& dT=bVec[_order-1];
    T sum(0);

    for(int i=1;i<_order;i++)
      {
	bVec[i-1]=(bVec[i-1]-(*this)(i,_order)*dT)/(*this)(i,i);
	sum+=(*this)(i,i)*bVec[i-1];
      }

    dT=(bVec[_order-1]-sum)/(*this)(_order,_order);

  }

  void DummySolve()
  {
    TArray bVec(_order);
    bVec=1.;
    //Replaces bVec with solution vector.
    
    T ani;
    T ain;
    T aii;
    T alpha;
    T beta;

    alpha=0.;
    beta=0.;
    
    for(int i=1;i<_order;i++)
      {
	ani=(*this)(_order,i);
	aii=(*this)(i,i);
	ain=(*this)(i,_order);
	alpha+=ani*bVec[i-1]/aii;
	beta+=ani*ain/aii;
      }

    T bn;
    bn=(bVec[_order-1]-alpha)/((*this)(_order,_order)-beta);
    
    for(int i=1;i<_order;i++)
      {
	ain=(*this)(i,_order);
	aii=(*this)(i,i);
	bVec[i-1]=(bVec[i-1]-ain*bn)/aii;
      }
    
    bVec[_order-1]=bn;
  }

  void SolveBotCol(TArray& bVec)
  {

    T sum(0);
    for(int i=1;i<_order+1;i++)
      sum+=getElement(_order,i);
    
    bVec[_order-1]/=sum;
    
  }
 
  void zero()
  {
    for(int i=0;i<_elements;i++)
      _values[i]=T_Scalar(0);
  }

  void ones()
  {
    for(int i=0;i<_elements;i++)
      _values[i]=T_Scalar(1);
  }

  T getTraceAbs()
  {
    T trace=0.;
    for(int i=0;i<_order;i++)
      trace+=fabs(_values[i]);
    return trace;
  }

  void multiply(const TArray& x, TArray& b)
  { //Ax=b
      int len=x.getLength();
      if(len==_order)
	{
	  b.zero();
	  for(int i=0;i<_order-1;i++)
	    b[i]=_values[i]*x[i]+_values[i-1+_order]*x[_order-1];
	  
	  b[_order-1]+=_values[_order-1]*x[_order-1];
	  for(int i=0;i<_order-1;i++)
	    b[_order-1]+=_values[2*_order+i-2]*x[i];
	}
      else
	throw CException("Array length does not match matrix order!");
    }
  
 private:
  int _elements;
  TArray _values;
  int _order;

};

#endif
