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
    
    for(int i=0;i<_SelfToOther.size();i++)
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

    for(int i=0;i<_SelfToSelf.size();i++)
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

    delete _empty.first;
    delete _empty.second;

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
  
 private:

  KSConnectivity(const KSConnectivity&);
  SelfToOther _SelfToOther;
  SelfToSelf _SelfToSelf;
  CouplingPair _empty;
  int _colLen;

};

#endif
