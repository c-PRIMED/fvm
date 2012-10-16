#ifndef _SCATTERINGKERNEL_H_
#define _SCATTERINGKERNEL_H_

#include "Array.h"
#include "KSConnectivity.h"
#include <iostream>
#include <fstream>

template <class T>
class ScatteringKernel
{
 public:

  ScatteringKernel() {}
  virtual void updateSource(const int cell)=0;

};

template <class T>
class GrapheneScatteringKernel : public ScatteringKernel<T>
{
  
  typedef Array<T> TArray;
  typedef Array<int> IntArray;

 public:
 GrapheneScatteringKernel(const int modeNum) :
  _type1Collisions(modeNum),
    _type2Collisions(modeNum)
    {
      _type1Collisions.setColumnLength();
      _type2Collisions.setColumnLength();
    }

  void ReadType1(const char* NamePhonon2, const char* NamePhonon3)
  {
    ifstream fp_in2;
    ifstream fp_in3;
    fp_in2.open(NamePhonon2,ifstream::in);
    fp_in3.open(NamePhonon3,ifstream::in);
    const int rowLen=_type1Collisions.getSelfSize();
    
    int row2;
    int nnz2;
    int row3;
    int nnz3;
    int index2;
    int index3;
    T dkl;
    T phi;

    for(int i=0;i<rowLen;i++)
      {
	fp_in2>>row2;
	fp_in2>>nnz2;
	fp_in3>>row3;
	fp_in3>>nnz3;

	if(nnz2!=0)
	  {
	    _type1Collisions.makeSelf(row2,nnz2);
	    _type1Collisions.makeOther(row3,nnz3);

	    for(int j=0;j<nnz2;j++)
	      {
		fp_in2>>index2;
		fp_in2>>dkl;
		fp_in3>>index3;
		fp_in3>>phi;
		_type1Collisions.setSelf(row2,j,index2,dkl);
		_type1Collisions.setOther(row3,j,index3,phi);
	      }

	  }

      }

  }

  void ReadType2(const char* NamePhonon2, const char* NamePhonon3)
  {
    ifstream fp_in2;
    ifstream fp_in3;
    fp_in2.open(NamePhonon2,ifstream::in);
    fp_in3.open(NamePhonon3,ifstream::in);
    const int rowLen=_type1Collisions.getSelfSize();
    
    int row2;
    int nnz2;
    int row3;
    int nnz3;
    int index2;
    int index3;
    T dkl;
    T phi;

    for(int i=0;i<rowLen;i++)
      {
	fp_in2>>row2;
	fp_in2>>nnz2;
	fp_in3>>row3;
	fp_in3>>nnz3;

	if(nnz2!=0)
	  {
	    _type2Collisions.makeSelf(row2,nnz2);
	    _type2Collisions.makeOther(row3,nnz3);

	    for(int j=0;j<nnz2;j++)
	      {
		fp_in2>>index2;
		fp_in2>>dkl;
		fp_in3>>index3;
		fp_in3>>phi;
		_type2Collisions.setSelf(row2,j,index2,dkl);
		_type2Collisions.setOther(row3,j,index3,phi);
	      }

	  }

      }

  }

  void updateSourceTerm(const TArray& e, const TArray& w, TArray& S)
  {

    S.zero();
    for(int row=0;row<e.getLength();row++)
      {
	const IntArray& t1p2=_type1Collisions.getSelfIndices(row);
	const IntArray& t1p3=_type1Collisions.getOtherIndices(row);
	const IntArray& t2p2=_type2Collisions.getSelfIndices(row);
	const IntArray& t2p3=_type2Collisions.getOtherIndices(row);

	const TArray& t1dk1=_type1Collisions.getSelfCoeffs(row);
	const TArray& t1phi=_type1Collisions.getOtherCoeffs(row);
	const TArray& t21dk1=_type2Collisions.getSelfCoeffs(row);
	const TArray& t21phi=_type2Collisions.getOtherCoeffs(row);

      }

  }

 private:
  GrapheneScatteringKernel(const GrapheneScatteringKernel&);
  KSConnectivity<T> _type1Collisions;
  KSConnectivity<T> _type2Collisions;

};

#endif
