// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KSCONNECTIVITY_H_
#define _KSCONNECTIVITY_H_

#include "Array.h"

template<class T>
class KSConnectivity
{
 public:

  typedef Array<int> IntArray;
  typedef Array<T> TArray;
  typedef pair<IntArray*,TArray*> CouplingPair;
  typedef vector<CouplingPair*> SelfToOther;
  typedef vector<CouplingPair*> SelfToSelf;
  
 KSConnectivity(const int selfLength, const int otherLength):
  _selfSite(selfLength),
    _otherSite(otherLength),
    _SelfToOtherCoeffs(0),
    _SelfToSelfCoeffs(0),
    _SelfToOtherConn(_selfSite,_otherSite),
    _SelfToSelfConn(_selfSite,_selfSite),
    _selfNNZ(0),
    _otherNNZ(0)
      {}

  void emptyConnections()
  {
    _SelfToOtherCoeffs.zero();
    _SelfToSelfCoeffs.zero();
  }

  void initSelfCount() {_SelfToSelfConn.initCount();}
  void initOtherCount() {_SelfToOtherConn.initCount();}
  void addCountSelf(const int index, const int count)
  {
    _SelfToSelfConn.addCount(index, count);
    _selfNNZ+=count;
  }
  void addCountOther(const int index, const int count)
  {
    _SelfToOtherConn.addCount(index, count);
    _otherNNZ+=count;
  }
  void finishCountSelf() 
  {
    _SelfToSelfConn.finishCount();
    _SelfToSelfCoeffs.resize(_selfNNZ);
    _SelfToSelfCoeffs.zero();
  }
  void finishCountOther() 
  {
    _SelfToOtherConn.finishCount();
    _SelfToOtherCoeffs.resize(_otherNNZ);
    _SelfToOtherCoeffs.zero();
  }
  void addSelf(const int i, const int j, const T val)
  {_SelfToSelfCoeffs[_SelfToSelfConn.add(i,j)]=val;}
  void addOther(const int i, const int j, const T val)
  {_SelfToOtherCoeffs[_SelfToOtherConn.add(i,j)]=val;}
  void finishAddSelf() {_SelfToSelfConn.finishAdd();}
  void finishAddOther() {_SelfToOtherConn.finishAdd();}
  int getSelfCount(const int i) {return _SelfToSelfConn.getCount(i);}
  int getOtherCount(const int i) {return _SelfToOtherConn.getCount(i);}
  const IntArray& getSelfRow() {return _SelfToSelfConn.getRow();}
  const IntArray& getOtherRow() {return _SelfToOtherConn.getRow();}
  const IntArray& getSelfCol() {return _SelfToSelfConn.getCol();}
  const IntArray& getOtherCol() {return _SelfToOtherConn.getCol();}
  const TArray& getSelfCoeffs() {return _SelfToSelfCoeffs;}
  const TArray& getOtherCoeffs() {return _SelfToOtherCoeffs;}
  TArray& getNonConstOtherCoeffs() {return _SelfToOtherCoeffs;}
  int getSelfNNZ() {return _selfNNZ;}
  int getOtherNNZ() {return _otherNNZ;}

  void multiplySelf(const TArray& x, TArray& b) const
  {//b=this*x
    const int Arows=_selfSite.getSelfCount();
    b.zero();
    if(Arows==x.getLength() && Arows==b.getLength())
      {
	const IntArray& row=_SelfToSelfConn.getRow();
	const IntArray& col=_SelfToSelfConn.getCol();
	for(int i=0;i<Arows;i++)
	  {
	    for(int pos=row[i];pos<row[i+1];pos++)
		b[i]+=x[col[pos]]*_SelfToSelfCoeffs[pos];
	  }
      }
    else
      throw CException("Matrix size does not agree with vectors!");
  }

  void multiplyOther(const TArray& x, TArray& b) const
  {//b=this*x
    const int Arows=_otherSite.getSelfCount();
    b.zero();
    if(Arows==x.getLength() && Arows==b.getLength())
      {
	for(int i=0;i<Arows;i++)
	  {
	    const IntArray& row=_SelfToOtherConn.getRow();
	    const IntArray& col=_SelfToOtherConn.getCol();
	    for(int pos=row[i];pos<row[i+1];pos++)
	      b[i]+=x[col[pos]]*_SelfToOtherCoeffs[pos];
	  }
      }
    else
      throw CException("Matrix size does not agree with vectors!");
  }
  
  int getSelfSize() {return _selfSite.getSelfCount();}
  int getOtherSize() {return _otherSite.getSelfCount();}

  void copyFrom(KSConnectivity& from)
  {
    //Copy SelfToSelf connections first
    
    initSelfCount();
    int Arows=from.getSelfSize();
    for(int i=0;i<Arows;i++)
      addCountSelf(i,from.getSelfCount(i));

    finishCountSelf();

    for(int i=0;i<Arows;i++)
      {
	const IntArray& row=from.getSelfRow();
	const IntArray& col=from.getSelfCol();
	const TArray& coeff=from.getSelfCoeffs();
	for(int pos=row[i];pos<row[i+1];pos++)
	  {
	    const T fromCoeff=coeff[pos];
	    const int j=col[pos];
	    addSelf(i,j,fromCoeff);
	  }
      }

    finishAddSelf();

    //Copy SelfToOther connections now

    initOtherCount();
    for(int i=0;i<Arows;i++)
      addCountOther(i,from.getOtherCount(i));

    finishCountOther();

    for(int i=0;i<Arows;i++)
      {
	const IntArray& row=from.getOtherRow();
	const IntArray& col=from.getOtherCol();
	const TArray& coeff=from.getOtherCoeffs();
	for(int pos=row[i];pos<row[i+1];pos++)
	  {
	    const T fromCoeff=coeff[pos];
	    const int j=col[pos];
	    addOther(i,j,fromCoeff);
	  }
      }

    finishAddOther();

  }

  void multiplySelf(const T x)
  {
    for(int i=0;i<_SelfToSelfCoeffs.getLength();i++)
      _SelfToSelfCoeffs[i]*=x;
  }

  void multiplyOther(const T x)
  {
    for(int i=0;i<_SelfToOtherCoeffs.getLength();i++)
      _SelfToOtherCoeffs[i]*=x;
  }

  void addToSelf(KSConnectivity& added)
  {
    const int selfSize=getSelfSize();
    if(selfSize==added.getSelfSize())
      {
	KSConnectivity<T> myCopy(getSelfSize(),getSelfSize());
	myCopy.copyFrom(*this);
	//start the counting
	initSelfCount();
	for(int i=0;i<selfSize;i++)
	  {
	    TArray myExpand(1);
	    myCopy.expandMySelfSelf(i, myExpand);

	    TArray addExpand(1);
	    added.expandMySelfSelf(i, addExpand);

	    myExpand+=addExpand;

	    int newSize(0);
	    for(int j=0;j<selfSize;j++)
	      {
		if(fabs(myExpand[j])>0.)
		  newSize++;
	      }
	    addCountSelf(i,newSize);
	  }
	finishCountSelf();
	    
	//put in values
	for(int i=0;i<selfSize;i++)
	  {
	    TArray myExpand(1);
	    myCopy.expandMySelfSelf(i, myExpand);

	    TArray addExpand(1);
	    added.expandMySelfSelf(i, addExpand);

	    myExpand+=addExpand;

	    for(int j=0;j<selfSize;j++)
	      {
		if(fabs(myExpand[j])>0.)
		    addSelf(i,j,myExpand[j]);
	      }
	  }
	finishAddSelf();

      }
    else
      throw CException("addToSelf: Rows not the same size!");
  }

  void addToOther(KSConnectivity& added)
  {
    const int selfSize=getSelfSize();
    if(selfSize==added.getSelfSize())
      {
	if(getOtherSize()==added.getOtherSize())
	  {
	    KSConnectivity<T> myCopy(getSelfSize(),getOtherSize());
	    myCopy.copyFrom(*this);
	    //start the counting
	    initOtherCount();
	    for(int i=0;i<selfSize;i++)
	      {
		TArray myExpand(1);
		myCopy.expandMySelfOther(i, myExpand);

		TArray addExpand(1);
		added.expandMySelfOther(i, addExpand);

		myExpand+=addExpand;

		int newSize(0);
		for(int j=0;j<getOtherSize();j++)
		  {
		    if(fabs(myExpand[j])>0.)
		      newSize++;
		  }
		addCountOther(i,newSize);
	      }
	    finishCountOther();
	    
	    //put in values
	    for(int i=0;i<selfSize;i++)
	      {
		TArray myExpand(1);
		myCopy.expandMySelfOther(i, myExpand);

		TArray addExpand(1);
		added.expandMySelfOther(i, addExpand);

		myExpand+=addExpand;

		for(int j=0;j<getOtherSize();j++)
		  {
		    if(fabs(myExpand[j])>0.)
			addOther(i,j,myExpand[j]);
		  }
	      }
	    finishAddOther();
	  }
	else
	  throw CException("addToOther: Columns not the same size!");
      }
    else
      throw CException("addToOther: Rows not the same size!");
  }

  void expandMySelfSelf(const int i, TArray& ExpCoeff)
  {
    const int mysize=getSelfSize();
    ExpCoeff.resize(mysize);
    ExpCoeff.zero();
    const IntArray& row=_SelfToSelfConn.getRow();
    const IntArray& col=_SelfToSelfConn.getCol();

    for(int pos=row[i];pos<row[i+1];pos++)
      {
	const int j=col[pos];
	ExpCoeff[j]=_SelfToSelfCoeffs[pos];
      }
    
  }

  void expandMySelfOther(const int i, TArray& ExpCoeff)
  {
    const int mysize=getOtherSize();
    ExpCoeff.resize(mysize);
    ExpCoeff.zero();
    const IntArray& row=_SelfToOtherConn.getRow();
    const IntArray& col=_SelfToOtherConn.getCol();

    for(int pos=row[i];pos<row[i+1];pos++)
      {
	const int j=col[pos];
	ExpCoeff[j]=_SelfToOtherCoeffs[pos];
      }
    
  }
  
 private:

  KSConnectivity(const KSConnectivity&);
  StorageSite _selfSite;
  StorageSite _otherSite;
  TArray _SelfToOtherCoeffs;
  TArray _SelfToSelfCoeffs;
  CRConnectivity _SelfToOtherConn;
  CRConnectivity _SelfToSelfConn;
  int _selfNNZ;
  int _otherNNZ;

};

#endif
