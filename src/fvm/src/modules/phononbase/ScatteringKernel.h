#ifndef _SCATTERINGKERNEL_H_
#define _SCATTERINGKERNEL_H_

#include "Array.h"
#include "KSConnectivity.h"
#include <iostream>
#include <fstream>


template <class T>
class ScatteringKernel
{
  
  typedef Array<T> TArray;
  typedef Array<int> IntArray;

 public:
 ScatteringKernel(const int modeNum) :
  _type1Collisions(modeNum, modeNum),
    _type2Collisions(modeNum, modeNum)
    {}

  void ReadType1(const char* NamePhonon2, const char* NamePhonon3)
  {

    cout<<"Reading type 1 collisions..."<<endl;
    ifstream fp_in2;
    ifstream fp_in3;
    fp_in2.open(NamePhonon2,ifstream::in | ifstream::binary);
    fp_in3.open(NamePhonon3,ifstream::in | ifstream::binary);
    const int rowLen=_type1Collisions.getSelfSize();
    
    int row2(0);
    int nnz2(0);
    int row3(0);
    int nnz3(0);
    int index2(0);
    int index3(0);
    T dkl(0);
    T phi(0);

    int* row2p(&row2);
    int* nnz2p(&nnz2);
    int* row3p(&row3);
    int* nnz3p(&nnz3);
    int* index2p(&index2);
    int* index3p(&index3);
    T* dklp(&dkl);
    T* phip(&phi);

    _type1Collisions.initSelfCount();
    _type1Collisions.initOtherCount();

    cout<<"Counting interactions..."<<endl;

    for(int i=0;i<rowLen;i++)
      {
	
	if(i % 100==0)
	  cout<<"Row: "<<i<<endl;

	fp_in2.read((char *)row2p,sizeof(int));
	fp_in2.read((char *)nnz2p,sizeof(int));
	fp_in3.read((char *)row3p,sizeof(int));
	fp_in3.read((char *)row3p,sizeof(int));

	for(int j=0;j<nnz2;j++)
	  {
	    fp_in2.read((char *)index2p,sizeof(int));
	    fp_in2.read((char *)dklp,sizeof(T));
	    fp_in3.read((char *)index3p,sizeof(int));
	    fp_in3.read((char *)phip,sizeof(T));
	    _type1Collisions.addCountSelf(i,1);
	    _type1Collisions.addCountOther(i,1);
	  }
      }

    _type1Collisions.finishCountSelf();
    _type1Collisions.finishCountOther();

    fp_in2.seekg(0, ios::beg);
    fp_in3.seekg(0, ios::beg);

    cout<<"Entering values..."<<endl;

    for(int i=0;i<rowLen;i++)
      {

	if(i % 100==0)
	  cout<<"Row: "<<i<<endl;

	fp_in2>>row2;
	fp_in2>>nnz2;
	fp_in3>>row3;
	fp_in3>>nnz3;

	for(int j=0;j<nnz2;j++)
	  {
	    fp_in2>>index2;
	    fp_in2>>dkl;
	    fp_in3>>index3;
	    fp_in3>>phi;
	    _type1Collisions.addSelf(i,index2,dkl);
	    _type1Collisions.addOther(i,index3,phi);
	  }
      }

    _type1Collisions.finishAddSelf();
    _type1Collisions.finishAddOther();

    cout<<"Type 1 complete."<<endl;

  }
  /*
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

  */

  void updateSourceTerm(const TArray& e, const TArray& w, TArray& S)
  {

    const T hbarJoule=1.054571726e-34;
    const T JouleToeV=6.24150974e18;
    const T Acell=5.378395621705545e-20;
    

    S.zero();
    const int Rows=S.getLength();
    const IntArray& t1p2row=_type1Collisions.getSelfRow();
    const IntArray& t1p3row=_type1Collisions.getOtherRow();
    const IntArray& t1p2col=_type1Collisions.getSelfCol();
    const IntArray& t1p3col=_type1Collisions.getOtherCol();
    const IntArray& t2p2row=_type2Collisions.getSelfRow();
    const IntArray& t2p3row=_type2Collisions.getOtherRow();
    const IntArray& t2p2col=_type2Collisions.getSelfCol();
    const IntArray& t2p3col=_type2Collisions.getOtherCol();

    const TArray& t1dkl=_type1Collisions.getSelfCoeffs();
    const TArray& t1phi=_type1Collisions.getOtherCoeffs();
    const TArray& t2dkl=_type2Collisions.getSelfCoeffs();
    const TArray& t2phi=_type2Collisions.getOtherCoeffs();

    //type 1 collisions
    for(int i=0;i<Rows;i++)
      {
	for(int pos=t1p2row[i];pos<t1p2row[i+1];pos++)
	  {
	    const int j2=t1p2col[pos];
	    const int j3=t1p3col[pos];
	    const T dkl=t1dkl[pos];
	    const T phi=t1phi[pos];
	    const T n1=e[i]/w[i];
	    const T n2=e[j2]/w[j2];
	    const T n3=e[j3]/w[j3];
	    S+=dkl*phi*(n1*(n3-n2)+n3*(n2+1));
	  }
      }

    //type 2 collisions
    for(int i=0;i<Rows;i++)
      {
	for(int pos=t2p2row[i];pos<t2p2row[i+1];pos++)
	  {
	    const int j2=t2p2col[pos];
	    const int j3=t2p3col[pos];
	    const T dkl=t2dkl[pos];
	    const T phi=t2phi[pos];
	    const T n1=e[i]/w[i];
	    const T n2=e[j2]/w[j2];
	    const T n3=e[j3]/w[j3];
	    S+=dkl*phi*(n1*(n3-n2)+n3*(n2+1));
	  }
      }

    const T preFac=Acell/6.283185307/hbarJoule*JouleToeV;
    for(int i=0;i<Rows;i++)
      S[i]*=preFac*w[i];

  }

 private:
  ScatteringKernel(const ScatteringKernel&);
  KSConnectivity<T> _type1Collisions;
  KSConnectivity<T> _type2Collisions;

};

#endif
