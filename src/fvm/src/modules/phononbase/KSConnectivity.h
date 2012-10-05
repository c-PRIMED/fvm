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
  
 KSConnectivity(const int selfLength):
  _SelfToOther(selfLength, NULL),
    _SelfToSelf(selfLength, NULL),
    _empty(),
    _colLen(-1)
      {
	IntArray* emptyInt=new IntArray(0);
	TArray* emptyT=new TArray(0);
	_empty.first=emptyInt;
	_empty.second=emptyT;
      }

  void emptyConnections()
  {
    
    for(int i=0;i<int(_SelfToOther.size());i++)
      {
	if(_SelfToOther[i]!=NULL)
	  {
	    IntArray* iptr=(_SelfToOther[i])->first;
	    delete iptr;
	    delete (_SelfToOther[i])->second;
	    delete _SelfToOther[i];
	  }
      }

    _SelfToOther.clear();

    for(int i=0;i<int(_SelfToSelf.size());i++)
      {
	if(_SelfToSelf[i]!=NULL)
	  {
	    IntArray* iptr=(_SelfToSelf[i])->first;
	    delete iptr;
	    delete (_SelfToSelf[i])->second;
	    delete _SelfToSelf[i];
	  }
      }
    
    _SelfToSelf.clear();

    //delete _empty.first;
    //delete _empty.second;

  }

  void setColumnLength(const int len) {_colLen=len;}

  void makeOther(const int index, const int length)
  {
    if(_colLen!=-1)
      {
	if(_SelfToOther[index]==NULL)
	  {
	    TArray* newOtherCoeffs=new TArray(length);
	    *newOtherCoeffs=-1.;
	    IntArray* newOtherIndices=new IntArray(length);
	    *newOtherIndices=-1;
	    CouplingPair* newOtherCoupling=new CouplingPair();
	    newOtherCoupling->first=newOtherIndices;
	    newOtherCoupling->second=newOtherCoeffs;
	    _SelfToOther[index]=newOtherCoupling;
	  }
	else
	  {
	    cout<<"Index "<<index<<endl;
	    throw CException("Other Connectivity Already Set!");
	  }
      }
    else
      throw CException("Have not set column length!");
  }

  void makeSelf(const int index, const int length)
  {
    if(_colLen!=-1)
      {
	if(_SelfToSelf[index]==NULL)
	  {
	    TArray* newSelfCoeffs=new TArray(length);
	    *newSelfCoeffs=-1.;
	    IntArray* newSelfIndices=new IntArray(length);
	    *newSelfIndices=-1;
	    CouplingPair* newSelfCoupling=new CouplingPair();
	    newSelfCoupling->first=newSelfIndices;
	    newSelfCoupling->second=newSelfCoeffs;
	    _SelfToSelf[index]=newSelfCoupling;
	  }
	else
	  {
	    cout<<"Index "<<index<<endl;
	    throw CException("Self Connectivity Already Set!");
	  }
      }
    else
      throw CException("Have not set column length!");
  }

  void setOther(const int SelfIndex, const int position, const int OtherIndex, const T coeff)
  {
    if(_SelfToOther[SelfIndex]!=NULL)
      {
	CouplingPair& Cpair=*_SelfToOther[SelfIndex];
	IntArray& positArray=(*Cpair.first);
	TArray& coeffArray=(*Cpair.second);
	const int length=positArray.getLength();
	if(position<length)
	  {
	    if(positArray[position]==-1)
	      {
		positArray[position]=OtherIndex;
		coeffArray[position]=coeff;
	      }
	    else
	      {
		cout<<"Index "<<SelfIndex<<"Position "<<position<<endl;
		throw CException("Other Connection has already been set!");
	      }
	  }
	else
	  {
	    cout<<"Index "<<SelfIndex<<"Position "<<position<<endl;
	    throw CException("Out of bounds for Other!");
	  }
      }
    else
      {
	cout<<"Index "<<SelfIndex<<endl;
	throw CException("Other Connectivity Not Made Yet!");
      }
  }

  void resetOther(const int SelfIndex, const int position, const int OtherIndex, const T coeff)
  {
    if(_SelfToOther[SelfIndex]!=NULL)
      {
	CouplingPair& Cpair=*_SelfToOther[SelfIndex];
	IntArray& positArray=(*Cpair.first);
	TArray& coeffArray=(*Cpair.second);
	const int length=positArray.getLength();
	if(position<length)
	  {
	    positArray[position]=OtherIndex;
	    coeffArray[position]=coeff;
	  }
	else
	  {
	    cout<<"Index "<<SelfIndex<<"Position "<<position<<endl;
	    throw CException("Out of bounds for Other!");
	  }
      }
    else
      {
	cout<<"Index "<<SelfIndex<<endl;
	throw CException("Other Connectivity Not Made Yet!");
      }
  }

  void setSelf(const int SelfIndex, const int position, const int OtherIndex, const T coeff)
  {
    if(_SelfToSelf[SelfIndex]!=NULL)
      {
	CouplingPair& Cpair=*_SelfToSelf[SelfIndex];
	IntArray& positArray=(*Cpair.first);
	TArray& coeffArray=(*Cpair.second);
	const int length=positArray.getLength();
	if(position<length)
	  {
	    if(positArray[position]==-1)
	      {
		positArray[position]=OtherIndex;
		coeffArray[position]=coeff;
	      }
	    else
	      {
		cout<<"Index "<<SelfIndex<<"Position "<<position<<endl;
		throw CException("Self Connection has already been set!");
	      }
	  }
	else
	  {
	    cout<<"Index "<<SelfIndex<<"Position "<<position<<endl;
	    throw CException("Out of bounds for Self!");
	  }
      }
    else
      {
	cout<<"Index "<<SelfIndex<<endl;
	throw CException("Self Connectivity Not Made Yet!");
      }
  }

  void resetSelf(const int SelfIndex, const int position, const int OtherIndex, const T coeff)
  {
    if(_SelfToSelf[SelfIndex]!=NULL)
      {
	CouplingPair& Cpair=*_SelfToSelf[SelfIndex];
	IntArray& positArray=(*Cpair.first);
	TArray& coeffArray=(*Cpair.second);
	const int length=positArray.getLength();
	if(position<length)
	  {
	    positArray[position]=OtherIndex;
	    coeffArray[position]=coeff;
	  }
	else
	  {
	    cout<<"Index "<<SelfIndex<<"Position "<<position<<endl;
	    throw CException("Out of bounds for Self!");
	  }
      }
    else
      {
	cout<<"Index "<<SelfIndex<<endl;
	throw CException("Self Connectivity Not Made Yet!");
      }
  }

  const IntArray& getOtherIndices(const int index) const
  {
    if(_SelfToOther[index]!=NULL)
      return *(_SelfToOther[index]->first);
    return *(_empty.first);
  }

  const IntArray& getSelfIndices(const int index) const
  {
    if (_SelfToSelf[index]!=NULL)
      return *(_SelfToSelf[index]->first);
    return *(_empty.first);
  }

  const TArray& getOtherCoeffs(const int index) const
  {
    if(_SelfToOther[index]!=NULL)
      return *(_SelfToOther[index]->second);
    return *(_empty.second);
  }

  const TArray& getSelfCoeffs(const int index) const
  {
    if(_SelfToSelf[index]!=NULL)
      return *(_SelfToSelf[index]->second);
    return *(_empty.second);
  }

  void multiplySelf(const TArray& x, TArray& b) const
  {//b=Ax
    const int Arows=_SelfToSelf.size();
    b.zero();
    if(Arows==x.getLength() && Arows==b.getLength())
      {
	for(int i=0;i<Arows;i++)
	  {
	    const IntArray& cIndex=getSelfIndices(i);
	    const TArray& cCoeff=getSelfCoeffs(i);
	    for(int j=0;j<cIndex.getLength();j++)
	      b[i]+=x[cIndex[j]]*cCoeff[j];
	  }
      }
    else
      throw CException("Matrix size does not agree with vectors!");
  }

  void multiplyOther(const TArray& x, TArray& b) const
  {//b=Ax
    const int Arows=_SelfToOther.size();
    b.zero();
    if(Arows==b.getLength())
      {
	for(int i=0;i<Arows;i++)
	  {
	    const IntArray& cIndex=getOtherIndices(i);
	    const TArray& cCoeff=getOtherCoeffs(i);
	    for(int j=0;j<cIndex.getLength();j++)
	      b[i]+=x[cIndex[j]]*cCoeff[j];
	  }
      }
    else
      throw CException("Matrix size does not agree with vectors!");
  }
  
  int getSelfSize() {return _SelfToSelf.size();}
  int getOtherSize() {return _SelfToOther.size();}
  bool isSelfNull(const int i) {return _SelfToSelf[i]==NULL;}
  bool isOtherNull(const int i) {return _SelfToOther[i]==NULL;}
  int getColSize() {return _colLen;}

  void copyFrom(KSConnectivity& from)
  {
    emptyConnections();

    //Copy SelfToSelf connections first
    
    int Arows=from.getSelfSize();
    for(int i=0;i<Arows;i++)
      {
	if(from.isSelfNull(i))
	  {
	    _SelfToSelf.push_back(NULL);
	  }
	else
	  {
	    const IntArray& cIndex=from.getSelfIndices(i);
	    const TArray& cCoeff=from.getSelfCoeffs(i);
	    const int newSize=cIndex.getLength();
	    IntArray* newIntsPtr=new IntArray(newSize);
	    TArray* newCoeffsPtr=new TArray(newSize);
	    *newIntsPtr=cIndex;
	    *newCoeffsPtr=cCoeff;
	    CouplingPair* newPairPtr=new CouplingPair(newIntsPtr, newCoeffsPtr);
	    _SelfToSelf.push_back(newPairPtr);
	  }
      }

    //Copy SelfToOther connections now

    Arows=from.getOtherSize();
    for(int i=0;i<Arows;i++)
      {
	if(from. isOtherNull(i))
	  {
	    _SelfToOther.push_back(NULL);
	  }
	else
	  {
	    const IntArray& cIndex=from.getOtherIndices(i);
	    const TArray& cCoeff=from.getOtherCoeffs(i);
	    const int newSize=cIndex.getLength();
	    IntArray* newIntsPtr=new IntArray(newSize);
	    TArray* newCoeffsPtr=new TArray(newSize);
	    *newIntsPtr=cIndex;
	    *newCoeffsPtr=cCoeff;
	    CouplingPair* newPairPtr=new CouplingPair(newIntsPtr, newCoeffsPtr);
	    _SelfToOther.push_back(newPairPtr);
	  }
      }

  }

  void multiplySelf(const T x)
  {
    const int Arows=_SelfToSelf.size();
    for(int i=0;i<Arows;i++)
      {
	if(!isSelfNull(i))
	  {
	    TArray& cCoeff=getSelfCoeffsPriv(i);
	    for(int j=0;j<cCoeff.getLength();j++)
	      cCoeff[j]*=x;
	  }
      }
  }

  void multiplyOther(const T x)
  {
    const int Arows=_SelfToOther.size();
    for(int i=0;i<Arows;i++)
      {
	if(!isOtherNull(i))
	  {
	    TArray& cCoeff=getOtherCoeffsPriv(i);
	    for(int j=0;j<cCoeff.getLength();j++)
	      cCoeff[j]*=x;
	  }
      }
  }

  void addToSelf(KSConnectivity& added)
  {
    const int selfSize=getSelfSize();
    if(selfSize==added.getSelfSize())
      {
	if(getColSize()==added.getColSize())
	  {
	    for(int i=0;i<selfSize;i++)
	      {
		if(isSelfNull(i) && !added.isSelfNull(i))
		  {
		    const IntArray& cIndex=added.getSelfIndices(i);
		    const TArray& cCoeff=added.getSelfCoeffs(i);
		    const int newSize=cIndex.getLength();
		    IntArray* newIntsPtr=new IntArray(newSize);
		    TArray* newCoeffsPtr=new TArray(newSize);
		    *newIntsPtr=cIndex;
		    *newCoeffsPtr=cCoeff;
		    CouplingPair* newPairPtr=new CouplingPair(newIntsPtr, newCoeffsPtr);
		    _SelfToSelf[i]=newPairPtr;
		  }
		else if(!isSelfNull(i) && !added.isSelfNull(i))
		  {
		    const IntArray& myIndex=getSelfIndices(i);
		    const TArray& myCoeff=getSelfCoeffs(i);
		    TArray myExpand(1);
		    expandArray(myIndex, myCoeff, myExpand);

		    const IntArray& addIndex=added.getSelfIndices(i);
		    const TArray& addCoeff=added.getSelfCoeffs(i);
		    TArray addExpand(1);
		    expandArray(addIndex, addCoeff, addExpand);

		    myExpand+=addExpand;

		    int newSize(0);
		    for(int j=0;j<_colLen;j++)
		      {
			if(fabs(myExpand[j])>0.)
			  newSize++;
		      }
		    
		    IntArray* newIntsPtr=new IntArray(newSize);
		    TArray* newCoeffsPtr=new TArray(newSize);

		    newSize=0;
		    for(int j=0;j<_colLen;j++)
		      {
			if(fabs(myExpand[j])>0.)
			  {
			    (*newIntsPtr)[newSize]=j;
			    (*newCoeffsPtr)[newSize]=myExpand[j];
			    newSize++;
			  }
		      }

		    CouplingPair* newPairPtr=new CouplingPair(newIntsPtr, newCoeffsPtr);
		    delete _SelfToSelf[i];
		    _SelfToSelf[i]=newPairPtr;
		  }
	      }
	  }
	else
	  throw CException("addToSelf: Columns not the same size!");
      }
    else
      throw CException("addToSelf: Rows not the same size!");
  }

  void addToOther(KSConnectivity& added)
  {
    const int otherSize=getOtherSize();
    if(otherSize==added.getOtherSize())
      {
	if(getColSize()==added.getColSize())
	  {
	    for(int i=0;i<otherSize;i++)
	      {
		if(isOtherNull(i) && !added.isOtherNull(i))
		  {
		    const IntArray& cIndex=added.getOtherIndices(i);
		    const TArray& cCoeff=added.getOtherCoeffs(i);
		    const int newSize=cIndex.getLength();
		    IntArray* newIntsPtr=new IntArray(newSize);
		    TArray* newCoeffsPtr=new TArray(newSize);
		    *newIntsPtr=cIndex;
		    *newCoeffsPtr=cCoeff;
		    CouplingPair* newPairPtr=new CouplingPair(newIntsPtr, newCoeffsPtr);
		    _SelfToOther[i]=newPairPtr;
		  }
		else if(!isOtherNull(i) && !added.isOtherNull(i))
		  {
		    const IntArray& myIndex=getOtherIndices(i);
		    const TArray& myCoeff=getOtherCoeffs(i);
		    TArray myExpand(1);
		    expandArray(myIndex, myCoeff, myExpand);

		    const IntArray& addIndex=added.getOtherIndices(i);
		    const TArray& addCoeff=added.getOtherCoeffs(i);
		    TArray addExpand(1);
		    expandArray(addIndex, addCoeff, addExpand);

		    myExpand+=addExpand;

		    int newSize(0);
		    for(int j=0;j<_colLen;j++)
		      {
			if(fabs(myExpand[j])>0.)
			  newSize++;
		      }
		    
		    IntArray* newIntsPtr=new IntArray(newSize);
		    TArray* newCoeffsPtr=new TArray(newSize);

		    newSize=0;
		    for(int j=0;j<_colLen;j++)
		      {
			if(fabs(myExpand[j])>0.)
			  {
			    (*newIntsPtr)[newSize]=j;
			    (*newCoeffsPtr)[newSize]=myExpand[j];
			    newSize++;
			  }
		      }

		    CouplingPair* newPairPtr=new CouplingPair(newIntsPtr, newCoeffsPtr);
		    delete _SelfToOther[i];
		    _SelfToOther[i]=newPairPtr;
		  }
	      }
	  }
	else
	  throw CException("addToOther: Columns not the same size!");
      }
    else
      throw CException("addToOther: Rows not the same size!");
  }

  void expandArray(const IntArray& indices, const TArray& compressed, TArray& expanded)
  {
    expanded.resize(_colLen);
    expanded.zero();
    const int len=indices.getLength();
    
    for(int i=0;i<len;i++)
      expanded[indices[i]]=compressed[i];
    
  }
  
 private:

  TArray& getSelfCoeffsPriv(const int index)
    {
      if(_SelfToSelf[index]!=NULL)
	return *(_SelfToSelf[index]->second);
      return *(_empty.second);
    }

  TArray& getOtherCoeffsPriv(const int index)
    {
      if(_SelfToOther[index]!=NULL)
	return *(_SelfToOther[index]->second);
      return *(_empty.second);
    }

  KSConnectivity(const KSConnectivity&);
  SelfToOther _SelfToOther;
  SelfToSelf _SelfToSelf;
  CouplingPair _empty;
  int _colLen;

};

#endif
