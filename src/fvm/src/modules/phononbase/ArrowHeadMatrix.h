#ifndef _ARROWHEADMATRIX_H_
#define _ARROWHEADMATRIX_H_

#include "MatrixJML.h"
#include "Array.h"
#include "ArrayBase.h"
#include "NumType.h"
#include <math.h>
#include <omp.h>

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

    #pragma omp parallel for default(shared) private(i,ani,aii,ain)
    {
    for(int i=1;i<_order;i++)
      {
	ani=getElement(_order,i);
	aii=getElement(i,i);
	ain=getElement(i,_order);

	alpha+=ani*bVec[i-1]/aii;
	beta+=ani*ain/aii;
      }
    }    

    T bn;
    bn=(bVec[_order-1]-alpha)/(getElement(_order,_order)-beta);
    
    #pragma omp parallel for default(shared) private(i,aii,ain)
    {
    for(int i=1;i<_order;i++)
      {
	ain=getElement(i,_order);
	aii=getElement(i,i);
	bVec[i-1]=(bVec[i-1]-ain*bn)/aii;
      }
    }

    bVec[_order-1]=bn;
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
    {
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
