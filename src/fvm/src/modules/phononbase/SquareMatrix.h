// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _SQUAREMATRIX_H_
#define _SQUAREMATRIX_H_

#include "Array.h"
#include "MatrixJML.h"
#include <math.h>

template<class T>
class SquareMatrix : public MatrixJML<T>
{
 public:
  typedef Array<T> TArray;
  typedef shared_ptr<TArray> TArrPtr;
  typedef Array<int> IntArray;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;

 SquareMatrix(const int N):
  _order(N),
    _elements(N*N),
    _sorted(false),
    _pivotRows(N),
    _maxVals(N),
    _values(_elements)
    {_values.zero();}
  
  T& getElement(const int i, const int j) {return _values[(i-1)*_order+j-1];}
  T& operator()(const int i, const int j) {return _values[(i-1)*_order+j-1];}
  void zero() {_values.zero();}
  
  void Solve(TArray& bVec)
  {//Gaussian Elimination w/ scaled partial pivoting
   //replaces bVec with the solution vector. 

    SquareMatrix<T> LU(_order);
    (*this).makeCopy(LU);
    TArray bCpy(_order);
    IntArray l(_order);
    TArray s(_order);

    //find max values in each row if not done yet
    if(!_sorted)
      {
	for(int i=1;i<_order+1;i++)
	  {
	    l[i-1]=i;
	    s[i-1]=fabs((*this)(i,1));
	    for(int j=2;j<_order+1;j++)
	      {
		if(s[i-1]<fabs((*this)(i,j)))
		  s[i-1]=fabs((*this)(i,j));
	      }
	    if(s[i-1]==0)
	      {
		cout<<"Row: "<<i<<endl;
		throw CException("Matrix has row of zeros!");
	      }
	  }
	_maxVals=s;
	//cout<<"Max vals:"<<endl;
	//for(int i=0;i<_order;i++)
	//cout<<s[i]<<endl;
      }
    
    //Forward sweep
    if(!_sorted)
      {
	for(int i=1;i<_order;i++)
	  {
	    T rmax=fabs((*this)(l[i-1],i)/s[l[i-1]-1]);
	    //  cout<<"rmax: "<<rmax<<endl;
	    int newMax=-1;
	    for(int j=i+1;j<_order+1;j++)
	      {
		T r=fabs((*this)(l[j-1],i)/s[l[j-1]-1]);
		if(r>rmax)
		  {
		    //cout<<"switching,i,j,r:"<<i<<","<<j<<","<<r<<endl;
		    //cout<<"factors "<<(*this)(l[j-1],i)<<" "<<s[l[j-1]-1]<<endl;
		    //cout<<"mapped j "<<l[j-1]<<endl;
		    rmax=r;
		    newMax=j;
		  }
	      }

	    if(newMax!=-1)
	      {
		int temp=l[i-1];
		l[i-1]=newMax;
		l[newMax-1]=temp;
	      }

	    for(int j=i+1;j<_order+1;j++)
	      {
		T factor=LU(l[j-1],i)/LU(l[i-1],i);
		LU(l[j-1],i)=factor;
		bVec[l[j-1]-1]-=factor*bVec[l[i-1]-1];
		for(int k=i+1;k<_order+1;k++)
		  {
		    LU(l[j-1],k)-=LU(l[i-1],k)*factor;
		    T test=LU(l[j-1],k);
		    if(isnan(test)||isinf(test))
		      {
			cout<<"Denom: "<<LU(l[i-1],i)<<endl;
			cout<<"Num: "<<LU(l[j-1],i)<<endl;
			cout<<"first: "<<LU(l[j-1],k)<<endl;
			cout<<"second: "<<LU(l[i-1],k)<<endl;
			cout<<"factor: "<<factor<<endl;
			cout<<"test: "<<test<<endl;
			throw CException("test is nan");
		      }
		  }
	      }
	  }
	_pivotRows=l;
	_sorted=true;
	//cout<<"Order"<<endl;
	//for(int i=0;i<_order;i++)
	//  cout<<_pivotRows[i]<<endl;
	//cout<<endl;
      }
    else
      {
	for(int i=1;i<_order;i++)
	  {
	    for(int j=i+1;j<_order+1;j++)
	      {
		T factor=LU(_pivotRows[j-1],i)/LU(_pivotRows[i-1],i);
		LU(_pivotRows[j-1],i)=factor;
		bVec[_pivotRows[j-1]-1]-=factor*bVec[_pivotRows[i-1]-1];
		for(int k=i+1;k<_order+1;k++)
		  {
		    LU(_pivotRows[j-1],k)=LU(_pivotRows[j-1],k)
		      -LU(_pivotRows[i-1],k)*factor;
		  }
	      }
	  }
      }

    //back solve
    bVec[_pivotRows[_order-1]-1]=
      bVec[_pivotRows[_order-1]-1]/LU(_pivotRows[_order-1],_order);

    T sum=0.;
    for(int i=_order-1;i>0;i--)
      {
	sum=0.;
	for(int j=i+1;j<_order+1;j++)
	  sum-=LU(_pivotRows[i-1],j)*bVec[_pivotRows[j-1]-1];
	bVec[_pivotRows[i-1]-1]+=sum;
	bVec[_pivotRows[i-1]-1]=bVec[_pivotRows[i-1]-1]/LU(_pivotRows[i-1],i);
      }
    
    //reorder
    bCpy=bVec;
    for(int i=0;i<_order;i++)
      bVec[i]=bCpy[_pivotRows[i]-1];

  }

  void makeCopy(SquareMatrix<T>& o)
  {
    if(o._order!=this->_order)
      throw CException("Cannot copy matrices of different sizes!");

    o._sorted=this->_sorted;
    o._pivotRows=this->_pivotRows;
    o._maxVals=this->_maxVals;
    o._values=this->_values;
  }

  void printMatrix()
  {
    for(int i=1;i<_order+1;i++)
      {
	for(int j=1;j<_order+1;j++)
	  cout<<(*this)(i,j)<<" ";
	cout<<endl;
      }
    cout<<endl;
  }

  T getTraceAbs()
  {
    T trace=0.;
    for(int i=1;i<_order+1;i++)
      trace+=fabs((*this)(i,i));
    return trace;
  }

  void multiply(const TArray& x, TArray& b)
    {
      int lenx=x.getLength();
      int lenb=b.getLength();
      if(lenx==_order && lenb==_order)
	{
	  b.zero();
	  for(int i=1;i<_order+1;i++)
	    for(int j=1;j<_order+1;j++)
	      b[i-1]+=(*this)(i,j)*x[j-1];
	}
      else
	throw CException("Array length does not match matrix order!");
    }

  void testSolve()
  {
    TArray x(_order);
    TArray b(_order);

    cout<<"Correct Solution:"<<endl;
    for(int i=0;i<_order;i++)
      {
	x[i]=rand()/1.e9;
	cout<<"x["<<i<<"]="<<x[i]<<endl;
      }
    
    multiply(x,b);
    
    cout<<"Before"<<endl;
    for(int i=0;i<_order;i++)
      cout<<"b["<<i<<"]="<<b[i]<<endl;

    Solve(b);

    cout<<"After"<<endl;
    for(int i=0;i<_order;i++)
      cout<<"b["<<i<<"]="<<b[i]<<endl;
    cout<<endl;

    cout<<"Error"<<endl;
    for(int i=0;i<_order;i++)
      cout<<"b["<<i<<"]="<<b[i]-x[i]<<endl;
    cout<<endl;
    

    throw CException("finished with test");

  }

 private:
  const int _order;
  const int _elements;
  bool _sorted;
  IntArray _pivotRows;
  TArray _maxVals;
  TArray _values;
  

};


#endif
