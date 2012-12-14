#ifndef _SCATTERINGKERNEL_H_
#define _SCATTERINGKERNEL_H_

#include "Array.h"
#include "KSConnectivity.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

template<class T> class Kspace;

template <class T>
class ScatteringKernel
{
  
  typedef Array<T> TArray;
  typedef Array<int> IntArray;
  typedef Array<bool> BArray;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef KSConnectivity<T> Tksconn;

 public:
 ScatteringKernel(Kspace<T>& kspace) :
  _kspace(kspace),
    _type1Collisions(kspace.gettotmodes(),kspace.gettotmodes()),
    _type2Collisions(kspace.gettotmodes(),kspace.gettotmodes()),
    _maxPhi(0),
    _maxDkl(0)
      {_kspace.setScattKernel(*this);}

  void ReadType1(const string& NamePhonon2, const string& NamePhonon3, const T tol)
  {

    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K)

    cout<<"Reading type I collisions..."<<endl;
    /*
    ifstream fp_in2;
    ifstream fp_in3;
    fp_in2.open(NamePhonon2,ios::in | ios::binary);
    fp_in3.open(NamePhonon3,ios::in | ios::binary);
    */

    FILE* fp_in2(NULL);
    FILE* fp_in3(NULL);

    if(fp_in3)
      fclose(fp_in3);
    if(fp_in2)
      fclose(fp_in2);

    int numAtt(10);
    while(fp_in2==NULL && numAtt>0)
      {
	fp_in2=fopen(NamePhonon2.c_str(),"rb");
	if(fp_in2==NULL)
	  cout<<"Error reading 2 "<<numAtt<<endl;
	numAtt--;
      }

    numAtt=10;
    while(fp_in3==NULL && numAtt>0)
      {
	fp_in3=fopen(NamePhonon3.c_str(),"rb");
	if(fp_in3==NULL)
	  cout<<"Error reading 3 "<<numAtt<<endl;
	numAtt--;
      }

    const int rowLen=_type1Collisions.getSelfSize();
    
    int row2(-1);
    int nnz2(-1);
    int row3(-1);
    int nnz3(-1);
    int index2(-1);
    int index3(-1);
    T dkl(-1);
    T phi(-1);

    _type1Collisions.initSelfCount();
    _type1Collisions.initOtherCount();

    TArray& w(_kspace.getFreqArray());

    cout<<"Counting I interactions..."<<endl;

    for(int i=0;i<rowLen;i++)
      {
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);
	const T n1=1./(exp(hbar*w[row2]/kb/300.)-1);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    //const T err=fabs((n1+1)*(n2+1)*n3-n1*n2*(n3+1))/n1;
	    const T err=fabs(w[row2]+w[index2]-w[index3])/w[row2];

	    const int m1=row2%6;
	    const int m2=index2%6;
	    const int m3=index3%6;

	    int zCnt(0);

	    if((m1==0)||(m1==3))
	      zCnt++;
	    if((m2==0)||(m2==3))
	      zCnt++;
	    if((m3==0)||(m3==3))
	      zCnt++;

	    //zCnt=0;
	    if(err<tol && zCnt!=1 && zCnt!=3 && dkl>0.)
	      {
		_type1Collisions.addCountSelf(row2,1);
		_type1Collisions.addCountOther(row3,1);
		if(phi>_maxPhi)
		  _maxPhi=phi;
		if(dkl>_maxDkl)
		  _maxDkl=dkl;
	      }

	  }
      }

    _type1Collisions.finishCountSelf();
    _type1Collisions.finishCountOther();
    
    rewind(fp_in2);
    rewind(fp_in3);

    cout<<"Entering I values..."<<endl;

    for(int i=0;i<rowLen;i++)
      {
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);
	const T n1=1./(exp(hbar*w[row2]/kb/300.)-1);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    //const T err=fabs((n1+1)*(n2+1)*n3-n1*n2*(n3+1))/n1;
	    const T err=fabs(w[row2]+w[index2]-w[index3])/w[row2];

	    const int m1=row2%6;
	    const int m2=index2%6;
	    const int m3=index3%6;

	    int zCnt(0);

	    if((m1==0)||(m1==3))
	      zCnt++;
	    if((m2==0)||(m2==3))
	      zCnt++;
	    if((m3==0)||(m3==3))
	      zCnt++;

	    //zCnt=0;
	    if(err<tol && zCnt!=1 && zCnt!=3 && dkl>0.)
		{
		  _type1Collisions.addSelf(row2,index2,dkl);
		  _type1Collisions.addOther(row3,index3,phi);
		}
	  }
      }

    _type1Collisions.finishAddSelf();
    _type1Collisions.finishAddOther();

    fclose(fp_in3);
    fclose(fp_in2);
    
    cout<<"Type I complete: "<<_type1Collisions.getSelfNNZ()<<endl;

  }

  void ReadType2(const string& NamePhonon2, const string& NamePhonon3, const T tol)
  {

    cout<<"Reading type II collisions..."<<endl;
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K)

    FILE* fp_in2(NULL);
    FILE* fp_in3(NULL);

    if(fp_in3)
      fclose(fp_in3);
    if(fp_in2)
      fclose(fp_in2);

    int numAtt(10);
    while(fp_in2==NULL && numAtt>0)
      {
	fp_in2=fopen(NamePhonon2.c_str(),"rb");
	if(fp_in2==NULL)
	  cout<<"Error reading 2 "<<numAtt<<endl;
	numAtt--;
      }

    numAtt=10;
    while(fp_in3==NULL && numAtt>0)
      {
	fp_in3=fopen(NamePhonon3.c_str(),"rb");
	if(fp_in3==NULL)
	  cout<<"Error reading 3 "<<numAtt<<endl;
	numAtt--;
      }

    TArray& w(_kspace.getFreqArray());
    const int rowLen=_type2Collisions.getSelfSize();
    
    int row2(-1);
    int nnz2(-1);
    int row3(-1);
    int nnz3(-1);
    int index2(-1);
    int index3(-1);
    T dkl(-1);
    T phi(-1);

    _type2Collisions.initSelfCount();
    _type2Collisions.initOtherCount();

    cout<<"Counting II interactions..."<<endl;

    for(int i=0;i<rowLen;i++)
      {
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);
	const T n1=1./(exp(hbar*w[row2]/kb/300.)-1);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    //const T err=fabs(n3*n2*(n1+1)-n1*(n2+1)*(n3+1))/n1;
	    const T err=fabs(w[row2]-w[index2]-w[index3])/w[index2];

	    const int m1=row2%6;
	    const int m2=index2%6;
	    const int m3=index3%6;

	    int zCnt(0);

	    if((m1==0)||(m1==3))
	      zCnt++;
	    if((m2==0)||(m2==3))
	      zCnt++;
	    if((m3==0)||(m3==3))
	      zCnt++;
	    
	    //zCnt=0;
	    if(zCnt!=1 && zCnt!=3 && dkl>0. && err<tol)
	      {
		_type2Collisions.addCountSelf(row2,1);
		_type2Collisions.addCountOther(row3,1);
		if(phi>_maxPhi)
		  _maxPhi=phi;
		if(dkl>_maxDkl)
		  _maxDkl=dkl;
	      }
	  }
      }

    _type2Collisions.finishCountSelf();
    _type2Collisions.finishCountOther();
    
    rewind(fp_in2);
    rewind(fp_in3);

    cout<<"Entering II values..."<<endl;

    for(int i=0;i<rowLen;i++)
      {
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);
	const T n1=1./(exp(hbar*w[row2]/kb/300.)-1);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    //const T err=fabs(n3*n2*(n1+1)-n1*(n2+1)*(n3+1))/n1;
	    const T err=fabs(w[row2]-w[index2]-w[index3])/w[index2];

	    const int m1=row2%6;
	    const int m2=index2%6;
	    const int m3=index3%6;

	    int zCnt(0);

	    if((m1==0)||(m1==3))
	      zCnt++;
	    if((m2==0)||(m2==3))
	      zCnt++;
	    if((m3==0)||(m3==3))
	      zCnt++;

	    //zCnt=0;
	    if(zCnt!=1 && zCnt!=3 && dkl>0. && err<tol)
	      {
		_type2Collisions.addSelf(row2,index2,dkl);
		_type2Collisions.addOther(row3,index3,phi);
	      }
	  }
      }

    _type2Collisions.finishAddSelf();
    _type2Collisions.finishAddOther();
    fclose(fp_in3);
    fclose(fp_in2);
    
    cout<<"Type II complete: "<<_type2Collisions.getSelfNNZ()<<endl;

    normalize();

  }

  void addFreqs()
  {
    const T hbarJoule=1.054571726e-34;
    TArray& w(_kspace.getFreqArray());
    const int Rows=w.getLength();

    const IntArray& t1p2row=_type1Collisions.getSelfRow();
    const IntArray& t1p3row=_type1Collisions.getOtherRow();
    const IntArray& t1p2col=_type1Collisions.getSelfCol();
    const IntArray& t1p3col=_type1Collisions.getOtherCol();
    const IntArray& t2p2row=_type2Collisions.getSelfRow();
    const IntArray& t2p3row=_type2Collisions.getOtherRow();
    const IntArray& t2p2col=_type2Collisions.getSelfCol();
    const IntArray& t2p3col=_type2Collisions.getOtherCol();

    TArray& t1phi=_type1Collisions.getNonConstOtherCoeffs();
    TArray& t2phi=_type2Collisions.getNonConstOtherCoeffs();

    //type 1 collisions
    
    for(int i=0;i<Rows;i++)
      {
	for(int pos=t1p2row[i];pos<t1p2row[i+1];pos++)
	  {
	    const int j2=t1p2col[pos];
	    const int j3=t1p3col[pos];
	    T& phi=t1phi[pos];
	    phi=phi/w[i]/w[j2]/w[j3];
	  }
      }

    //type 2 collisions
    for(int i=0;i<Rows;i++)
      {
	for(int pos=t2p2row[i];pos<t2p2row[i+1];pos++)
	  {
	    const int j2=t2p2col[pos];
	    const int j3=t2p3col[pos];
	    T& phi=t2phi[pos];
	    const T w3=w[i]-w[j2];
	    phi=phi/w[i]/w[j2]/w[j3];
	  }
      }
    
  }
 
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
	    S[i]+=dkl*phi*(n1*(n3-n2)+n3*(n2+1));
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
	    S[i]+=dkl*phi*(n1*(n3-n2)+n3*(n2+1));
	  }
      }

    const T preFac=Acell/6.283185307/hbarJoule*JouleToeV;
    for(int i=0;i<Rows;i++)
      S[i]*=preFac*w[i];

  }

  void updateSourceTermTest(const T Tl)
  {

    const T hbarJoule=1.054571726e-34;
    const T hbar=6.582119e-16;
    const T JouleToeV=6.24150974e18;
    const T Acell=5.378395621705545e-20;
    const T cv=1e-9*1e-9;
    T temp=Tl;
    
    TArray S(_kspace.gettotmodes());
    S.zero();
    TArray e(_kspace.gettotmodes());
    e.zero();
    _kspace.getEquilibriumArray(e,Tl);
    TArray& w(_kspace.getFreqArray());

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
	    const T phi=t1phi[pos]*Acell;
	    const T n1=e[i]/hbar/w[i];
	    const T n2=e[j2]/hbar/w[j2];
	    const T n3=e[j3]/hbar/w[j3];
	    S[i]+=dkl*phi*((n1+1)*(n2+1)*n3-n1*n2*(n3+1));
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
	    const T phi=t2phi[pos]*Acell;
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    S[i]+=0.5*dkl*phi*(n3*n2*(n1+1)-n1*(n2+1)*(n3+1));
	  }
      }

    const T preFac=1./2./3.141592653/pow(hbarJoule,2)*_maxDkl*_maxPhi*cv;
    for(int i=0;i<Rows;i++)
      S[i]*=preFac*w[i]*hbar;

    T defect(0);
    for(int i=0;i<Rows;i++)
      defect+=S[i];

    for(int i=0;i<Rows;i++)
      e[i]+=S[i];

    int cnt(0);
    T esum(0);
    const int klen=_kspace.getlength();
    for(int k=0;k<klen;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    esum+=e[cnt]*dk3;
	    cnt++;
	  }
      }

    _kspace.calcTemp(temp,esum);
    cout<<"New Temp: "<<temp<<endl;

  }

  void updateSource(const int c, TArray& S)
  {

    const T hbarJoule=1.054571726e-34;
    const T hbar=6.582119e-16;
    const T Acell=5.378395621705545e-20;
    const T kb=8.617343e-5;  // (eV/K)
    T cOt(0);
    S.zero();
    const int Rows=S.getLength();
    TArray e(Rows);
    //TArray e0(Rows);
    _kspace.geteCellVals(c,e);
    //_kspace.gete0CellVals(c,e0);
    const TArray& w(_kspace.getFreqArray());
    const T Tlat(300.05);//Tl(_kspace.calcLatTemp(c));
    const T Tl(_kspace.calcLatTemp(c));
    TArray e0(Rows);
    //_kspace.gete0CellVals(c,e0);
    _kspace.getEquilibriumArray(e0,Tl);

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
	    const T w3=w[i]+w[j2];
	    const T w1=w[j3]-w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    const T n10=e0[i]/w[i]/hbar;
	    const T n20=e0[j2]/w[j2]/hbar;
	    const T n30=e0[j3]/w[j3]/hbar;
	    const T n012=1./(exp(hbar*w3/kb/Tl)-1.);
	    const T dn1=n1-n10;
	    const T dn2=n2-n20;
	    const T dn3=(n3-n30);//*n012/n30;
	    T Large(n10);
	    T e032(0);
	    T frac(0);
	    //const T frac=w[j3]*(1-w[j2]/w[j3])/w[i];
	    
	    if(w1>0.)
	      {
		e032=hbar*w1/(exp(hbar*w1/kb/Tl)-1.);
		frac=e032/e0[i];
	      }
	     
	    if(n20>Large)
	      Large=n20;
	    if(n30>Large)
	      Large=n30;

	    const T nsum=Large*(-dn1*dn2/Large+dn3/Large+dn1*dn3/Large
				+dn2*dn3/Large-dn2*n10/Large+dn3*n10/Large-dn1*n20/Large
				+dn3*n20/Large+dn1*n30/Large+dn2*n30/Large);
	    /*
	    if(j2>i)
	      {
		S[i]+=dkl*phi*nsum;
		S[j2]+=dkl*phi*nsum;
		S[j3]-=0.5*dkl*phi*nsum;
		//S[i]+=dkl*phi*((n1+1)*(n2+1)*n3-n1*n2*(n3+1));
	      }*/

	    S[i]+=dkl*phi*nsum;
	  }
      }
    
    T maxS(0);

    //type 2 collisions
    
    for(int i=0;i<Rows;i++)
      {
	for(int pos=t2p2row[i];pos<t2p2row[i+1];pos++)
	  {
	    const int j2=t2p2col[pos];
	    const int j3=t2p3col[pos];
	    const T dkl=t2dkl[pos];
	    const T phi=t2phi[pos];
	    const T w3=w[i]+w[j2];
	    const T w1=w[j3]-w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    const T n10=e0[i]/w[i]/hbar;
	    const T n20=e0[j2]/w[j2]/hbar;
	    const T n30=e0[j3]/w[j3]/hbar;
	    const T n012=1./(exp(hbar*w3/kb/Tl)-1.);
	    const T dn1=n1-n10;
	    const T dn2=n2-n20;
	    const T dn3=(n3-n30);//*n012/n30;
	    //const T frac=w[j3]*(1-w[j2]/w[j3])/w[i];
	    T e032(0);
	    T frac(0);
	    
	    if(w1>0.)
	      {
		e032=hbar*w1/(exp(hbar*w1/kb/Tl)-1.);
		frac=e032/e0[i];
	      }

	    T Large(n10);

	    if(n20>Large)
	      Large=n20;
	    if(n30>Large)
	      Large=n30;

	    const T nsum=Large*(-dn1/Large-dn1*dn2/Large-dn1*dn3/Large
				+dn2*dn3/Large-dn2*n10/Large-dn3*n10/Large
				-dn1*n20/Large+dn3*n20/Large-dn1*n30/Large
				+dn2*n30/Large);

	    S[i]+=0.5*dkl*phi*nsum;
	    //S[i]+=0.5*dkl*phi*(n3*n2*(n1+1)-n1*(n2+1)*(n3+1));
	  }
      }

    const T preFac=Acell/16./3.141592653*hbarJoule*_maxDkl*_maxPhi;
    const int klen=_kspace.getlength();
    for(int i=0;i<Rows;i++)
      {
	S[i]*=(preFac*w[i]*hbar);
	if(fabs(S[i])>fabs(maxS))
	   maxS=fabs(S[i]);
      }
    
    const T DK3=_kspace.getDK3();
    T defect(0);
    int cnt(0);
    T eq(0);
    const int cStart(_kspace.getGlobalIndex(c,0));
    int ind(cStart);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    defect+=((S[cnt]/maxS)*(dk3/DK3));
	    Tmode& mode=kv.getmode(m);
	    const T tau=mode.gettau();
	    cOt+=mode.calce0(Tl)*(dk3/DK3);
	    //const T e0=mode.calce0(Tl);
	    //eq+=(e0-e[ind])/tau*(dk3/DK3);
	    ind++;
	    cnt++;
	  }
      }

    cOt*=DK3;
    T f(maxS*defect/eq);
    defect*=maxS*DK3;
    //cout<<"Defect1: "<<defect<<endl;
    
    if(fabs(defect)>0.)
      {
	cnt=0;
	T d(defect);
	defect=0.;
	ind=cStart;
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		const T tau=mode.gettau();
		const T cp=mode.calcde0dT(Tl);
		const T e0(mode.calce0(Tl));
		//S[cnt]+=f*(e0-e[ind])/tau;
		S[cnt]=d*(S[cnt]/d-e0/cOt);
		defect+=((S[cnt]/maxS)*(dk3/DK3));
		ind++;
		cnt++;
	      }
	  }
	defect*=DK3*maxS;
	//cout<<"Defect2: "<<defect<<endl;
      }
    else
      {
	S.zero();
      }
    
  }

  void getTypeIsource(const int c, TArray& S, TArray& dS, const bool correct)
  {

    const T hbarJoule=1.054571726e-34;
    const T hbar=6.582119e-16;
    const T Acell=5.378395621705545e-20;
    const T kb=8.617343e-5;  // (eV/K)
    T cOt(0);
    S.zero();
    dS.zero();
    const int Rows=S.getLength();
    TArray Sl(Rows);
    TArray dSl(Rows);
    Sl.zero();
    dSl.zero();
    TArray e(Rows);
    //TArray e0(Rows);
    _kspace.geteCellVals(c,e);
    //_kspace.gete0CellVals(c,e0);
    const TArray& w(_kspace.getFreqArray());
    const T Tlat(300.0);//Tl(_kspace.calcLatTemp(c));
    const T Tl(_kspace.calcLatTemp(c));
    TArray e0(Rows);
    //_kspace.gete0CellVals(c,e0);
    _kspace.getEquilibriumArray(e0,Tl);
    const T DK3=_kspace.getDK3();

    const IntArray& t1p2row=_type1Collisions.getSelfRow();
    const IntArray& t1p3row=_type1Collisions.getOtherRow();
    const IntArray& t1p2col=_type1Collisions.getSelfCol();
    const IntArray& t1p3col=_type1Collisions.getOtherCol();
    const TArray& t1dkl=_type1Collisions.getSelfCoeffs();
    const TArray& t1phi=_type1Collisions.getOtherCoeffs();

    //type 1 collisions
    
    for(int i=0;i<Rows;i++)
      {
	const int k1=floor(i/6);
	const T dk31=_kspace.getkvol(k1).getdk3();
	for(int pos=t1p2row[i];pos<t1p2row[i+1];pos++)
	  {
	    const int j2=t1p2col[pos];
	    const int k2=floor(j2/6);
	    const T dk32=_kspace.getkvol(k2).getdk3();
	    const int j3=t1p3col[pos];
	    const int k3=floor(j3/6);
	    const T dk33=_kspace.getkvol(k3).getdk3();
	    const T dkl=t1dkl[pos];
	    const T phi=t1phi[pos];
	    const T w3=w[i]+w[j2];
	    const T w1=w[j3]-w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    const T n10=e0[i]/w[i]/hbar;
	    const T n20=e0[j2]/w[j2]/hbar;
	    const T n30=e0[j3]/w[j3]/hbar;
	    const T n012=1./(exp(hbar*w3/kb/Tl)-1.);
	    const T dn1=n1-n10;
	    const T dn2=n2-n20;
	    const T dn3=(n3-n30);//*n012/n30;
	    T Large(n10);
	    T e032(0);
	    T frac(0);
	    //const T frac=w[j3]*(1-w[j2]/w[j3])/w[i];

	    if(n20>Large)
	      Large=n20;
	    if(n30>Large)
	      Large=n30;

	    const T nsum=(-dn1*dn2+dn3+dn1*dn3
				+dn2*dn3-dn2*n10+dn3*n10-dn1*n20
				+dn3*n20+dn1*n30+dn2*n30);
	    
	    const T nsum1=(-dn1*n20+dn1*n30);
	    //const T nsum1=(n1+1)*(n20+1)*n30-n1*n20*(n30+1);
	    
	      
	    Sl[i]+=dkl*phi*nsum*w[i];
	    Sl[j2]+=dkl*phi*nsum*w[j2]*dk31/dk32;
	    Sl[j3]-=dkl*phi*nsum*w3*dk31/dk33;
	    //S[i]+=dkl*phi*((n1+1)*(n2+1)*n3-n1*n2*(n3+1));
	      
	  }
      }
    
    T maxS(0);
    const T preFac=Acell/16./3.141592653*hbarJoule*_maxPhi*hbar*_maxDkl;
    const int klen=_kspace.getlength();
    for(int i=0;i<Rows;i++)
      {
	Sl[i]*=preFac;//*w[i]*hbar);
	dSl[i]*=preFac;
	if(fabs(Sl[i])>fabs(maxS))
	   maxS=fabs(Sl[i]);
      }
    
    long double defect(0);
    int cnt(0);
    T eq(0);
    const int cStart(_kspace.getGlobalIndex(c,0));
    int ind(cStart);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    defect+=((Sl[cnt]/maxS)*(dk3/DK3));
	    //if(fabs(Sl[cnt])*(dk3/DK3)>100)
	    // Sl[cnt]=0.;
	    Tmode& mode=kv.getmode(m);
	    const T tau=mode.gettau();
	    cOt+=mode.calce0(Tl)*(dk3/DK3);
	    //const T e0=mode.calce0(Tl);
	    //eq+=(e0-e[ind])/tau*(dk3/DK3);
	    ind++;
	    cnt++;
	  }
      }

    cOt*=DK3;
    T f(maxS*defect/eq);
    defect*=maxS*DK3;
    //cout<<"Defect1: "<<defect<<endl;
    //cout<<"Temperature: "<<Tl<<endl;
    
    
    if(correct)
      {
	cnt=0;
	T d(defect);
	defect=0.;
	ind=cStart;
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		const T tau=mode.gettau();
		const T cp=mode.calcde0dT(Tl);
		const T e0(mode.calce0(Tl));
		//S[cnt]+=f*(e0-e[ind])/tau;
		//Sl[cnt]=d*(Sl[cnt]/d-e0/cOt);
		defect+=((Sl[cnt]/maxS)*(dk3/DK3));
		ind++;
		cnt++;
	      }
	  }
	defect*=DK3*maxS;
	//cout<<"Defect2: "<<defect<<endl;
      }

    for(int i=0;i<Rows;i++)
      {
	S[i]=Sl[i];
	//dS[i]=dSl[i];
      }

  }

  void updateSource2(const int c, TArray& S, TArray& dS)
  {

    const T hbarJoule=1.054571726e-34;
    const T hbar=6.582119e-16;
    const T Acell=5.378395621705545e-20;
    const T kb=8.617343e-5;  // (eV/K)
    T cOt(0);
    S.zero();
    dS.zero();
    const int Rows=S.getLength();
    TArray Sl(Rows);
    TArray dSl(Rows);
    Sl.zero();
    dSl.zero();
    TArray e(Rows);
    //TArray e0(Rows);
    _kspace.geteCellVals(c,e);
    //_kspace.gete0CellVals(c,e0);
    const TArray& w(_kspace.getFreqArray());
    const T Tlat(300.0);//Tl(_kspace.calcLatTemp(c));
    const T Tl(_kspace.calcLatTemp(c));
    TArray e0(Rows);
    _kspace.gete0CellVals(c,e0);
    //_kspace.getEquilibriumArray(e0,Tl);
    const T DK3=_kspace.getDK3();

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
	const int k1=floor(i/6);
	const T dk31=_kspace.getkvol(k1).getdk3();
	for(int pos=t1p2row[i];pos<t1p2row[i+1];pos++)
	  {
	    const int j2=t1p2col[pos];
	    const int k2=floor(j2/6);
	    const T dk32=_kspace.getkvol(k2).getdk3();
	    const int j3=t1p3col[pos];
	    const int k3=floor(j3/6);
	    const T dk33=_kspace.getkvol(k3).getdk3();
	    const T dkl=t1dkl[pos];
	    const T phi=t1phi[pos];
	    const T w3=w[i]+w[j2];
	    const T w1=w[j3]-w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    const T n10=e0[i]/w[i]/hbar;
	    const T n20=e0[j2]/w[j2]/hbar;
	    const T n30=e0[j3]/w[j3]/hbar;
	    const T n012=1./(exp(hbar*w3/kb/Tl)-1.);
	    const T dn1=n1-n10;
	    const T dn2=n2-n20;
	    const T dn3=n3-n30;
	    T e032(0);
	    T frac(0);
	    //const T frac=w[j3]*(1-w[j2]/w[j3])/w[i];

	    const T nsum=(-dn1*dn2+dn3+dn1*dn3
			  +dn2*dn3-dn2*n10+dn3*n10-dn1*n20
			  +dn3*n20+dn1*n30+dn2*n30);
	    
	    const T nsum1=(-dn1*n20+dn1*n30);
	    
	    const T ds1=-dn2+dn3-n20+n30;
	    const T ds2=-dn1+dn3-n10+n30;
	    const T ds3=1+dn1+dn2+n10+n20;
	    
	      
	    Sl[i]+=dkl*phi*nsum*w[i];//*dk33/dk31;
	    //Sl[j2]+=dkl*phi*nsum*w[j2]*dk33/dk32;
	    //Sl[j3]-=dkl*phi*nsum*w[j3]*dk31/dk33;
	    //Sl[j3]-=(dkl*phi*nsum*w[i]+dkl*phi*nsum*w[j2]);

	    dSl[i]+=dkl*phi*ds1*w[i];
	    dSl[j2]+=dkl*phi*ds2*w[j2]*dk31/dk32;
	    dSl[j3]-=dkl*phi*ds3*w[j3]*dk31/dk33;

	  }
      }
    
    T maxS(0);

    //type 2 collisions
    
    for(int i=0;i<Rows;i++)
      {
	const int k1=floor(i/6);
	const T dk31=_kspace.getkvol(k1).getdk3();
	for(int pos=t2p2row[i];pos<t2p2row[i+1];pos++)
	  {
	    const int j2=t2p2col[pos];
	    const int k2=floor(j2/6);
	    const T dk32=_kspace.getkvol(k2).getdk3();
	    const int j3=t2p3col[pos];
	    const int k3=floor(j3/6);
	    const T dk33=_kspace.getkvol(k3).getdk3();
	    const T dkl=t2dkl[pos];
	    const T phi=t2phi[pos];
	    const T w3=w[i]-w[j2];
	    const T w2=w[j3]-w[i];
	    const T w1=w[j3]-w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    const T n10=e0[i]/w[i]/hbar;
	    const T n20=e0[j2]/w[j2]/hbar;
	    const T n30=e0[j3]/w[j3]/hbar;
	    const T n012=1./(exp(hbar*w3/kb/Tl)-1.);
	    const T e012=hbar*w3/(exp(hbar*w3/kb/Tl)-1.);
	    const T dn1=n1-n10;
	    const T dn2=n2-n20;
	    const T dn3=n3-n30;
	    //const T frac=w[j3]*(1-w[j2]/w[j3])/w[i];

	    const T nsum=(-dn1-dn1*dn2-dn1*dn3
				+dn2*dn3-dn2*n10-dn3*n10
				-dn1*n20+dn3*n20-dn1*n30
				+dn2*n30);

	    const T ds1=-(1+dn2+dn3+n20+n30);
	    const T ds2=-dn1+dn3-n10+n30;
	    const T ds3=-dn1+dn2-n10+n20;

	    
	    const T nsum1=(-dn1-dn1*n20-dn1*n30);
	    const T nsum2=(-dn2*n10+dn2*n30);
	    const T nsum3=(-dn3*n10+dn3*n20);
	    

	    Sl[i]+=0.5*dkl*phi*nsum*w[i];
	    //Sl[i]+=0.5*dkl*phi*nsum*w[j2]+0.5*dkl*phi*nsum*w[j3];
	    //dSl[i]+=0.5*dkl*phi*ds1*w[i];
	    dSl[i]+=0.5*dkl*phi*ds1*w[j2]+0.5*dkl*phi*ds1*w[j3];

	    //Sl[j2]-=0.5*dkl*phi*nsum*w[j2]*dk31/dk32;
	    dSl[j2]-=0.5*dkl*phi*ds2*w[j2]*dk31/dk32;

	    //Sl[j3]-=0.5*dkl*phi*nsum*w[j3]*dk31/dk33;
	    dSl[j3]-=0.5*dkl*phi*ds3*w[j3]*dk31/dk33;

	  }
      }

    const T preFac=Acell/16./3.141592653*hbarJoule*_maxPhi*hbar*_maxDkl;
    const int klen=_kspace.getlength();
    for(int i=0;i<Rows;i++)
      {
	Sl[i]*=preFac;//*w[i]*hbar);
	dSl[i]*=preFac;
	if(fabs(Sl[i])>fabs(maxS))
	   maxS=fabs(Sl[i]);
      }
    
    long double defect(0);
    int cnt(0);
    T eq(0);
    const int cStart(_kspace.getGlobalIndex(c,0));
    int ind(cStart);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    defect+=((Sl[cnt]/maxS)*(dk3/DK3));
	    //if(fabs(Sl[cnt])*(dk3/DK3)>100)
	    // Sl[cnt]=0.;
	    Tmode& mode=kv.getmode(m);
	    const T tau=mode.gettau();
	    cOt+=mode.calce0(Tl)*(dk3/DK3);
	    //const T e0=mode.calce0(Tl);
	    //eq+=(e0-e[ind])/tau*(dk3/DK3);
	    ind++;
	    cnt++;
	  }
      }

    cOt*=DK3;
    T f(maxS*defect/eq);
    //cout<<"Defect1: "<<defect<<endl;
    //cout<<"Temperature: "<<Tl<<endl;

    if(fabs(defect)>0.)
      {
	defect*=maxS*DK3;
	cnt=0;
	T d(defect);
	defect=0.;
	ind=cStart;
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		const T tau=mode.gettau();
		const T cp=mode.calcde0dT(Tl);
		const T e0(mode.calce0(Tl));
		//S[cnt]+=f*(e0-e[ind])/tau;
		Sl[cnt]=d*(Sl[cnt]/d-e0/cOt);
		defect+=((Sl[cnt]/maxS)*(dk3/DK3));
		ind++;
		cnt++;
	      }
	  }
	//defect*=DK3*maxS;
	//cout<<"Defect2: "<<defect<<endl;
      }

    for(int i=0;i<Rows;i++)
      {
	S[i]=Sl[i];
	//dS[i]=dSl[i];
      }

  }

  void getTypeIIsource(const int c, TArray& S, TArray& dS, const bool correct)
  {

    const T hbarJoule=1.054571726e-34;
    const T hbar=6.582119e-16;
    const T Acell=5.378395621705545e-20;
    const T kb=8.617343e-5;  // (eV/K)
    T cOt(0);
    S.zero();
    dS.zero();
    const int Rows=S.getLength();
    TArray Sl(Rows);
    TArray dSl(Rows);
    Sl.zero();
    dSl.zero();
    TArray e(Rows);
    //TArray e0(Rows);
    _kspace.geteCellVals(c,e);
    //_kspace.gete0CellVals(c,e0);
    const TArray& w(_kspace.getFreqArray());
    const T Tlat(300.0);//Tl(_kspace.calcLatTemp(c));
    const T Tl(_kspace.calcLatTemp(c));
    TArray e0(Rows);
    //_kspace.gete0CellVals(c,e0);
    _kspace.getEquilibriumArray(e0,Tl);
    const T DK3=_kspace.getDK3();

    const IntArray& t2p2row=_type2Collisions.getSelfRow();
    const IntArray& t2p3row=_type2Collisions.getOtherRow();
    const IntArray& t2p2col=_type2Collisions.getSelfCol();
    const IntArray& t2p3col=_type2Collisions.getOtherCol();

    const TArray& t2dkl=_type2Collisions.getSelfCoeffs();
    const TArray& t2phi=_type2Collisions.getOtherCoeffs();
    
    T maxS(0);

    //type 2 collisions
    
    for(int i=0;i<Rows;i++)
      {
	const int k1=floor(i/6);
	const T dk31=_kspace.getkvol(k1).getdk3();
	for(int pos=t2p2row[i];pos<t2p2row[i+1];pos++)
	  {
	    const int j2=t2p2col[pos];
	    const int k2=floor(j2/6);
	    const T dk32=_kspace.getkvol(k2).getdk3();
	    const int j3=t2p3col[pos];
	    const int k3=floor(j3/6);
	    const T dk33=_kspace.getkvol(k3).getdk3();
	    const T dkl=t2dkl[pos];
	    const T phi=t2phi[pos];
	    const T w3=w[i]-w[j2];
	    const T w2=w[j3]-w[i];
	    const T w1=w[j3]-w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w[j3]/hbar;
	    const T n10=e0[i]/w[i]/hbar;
	    const T n20=e0[j2]/w[j2]/hbar;
	    const T n30=e0[j3]/w[j3]/hbar;
	    const T n012=1./(exp(hbar*w3/kb/Tl)-1.);
	    const T e012=hbar*w3/(exp(hbar*w3/kb/Tl)-1.);
	    const T dn1=n1-n10;
	    const T dn2=n2-n20;
	    const T dn3=(n3-n30);//*n012/n30;
	    //const T frac=w[j3]*(1-w[j2]/w[j3])/w[i];

	    T Large(n10);

	    if(n20>Large)
	      Large=n20;
	    if(n30>Large)
	      Large=n30;
	    
	    
	    const T nsum=(-dn1-dn1*dn2-dn1*dn3
				+dn2*dn3-dn2*n10-dn3*n10
				-dn1*n20+dn3*n20-dn1*n30
				+dn2*n30);

	    const T ds1=-(1+dn2+dn3-n20-n30);
	    const T ds2=-dn1+dn3-n10+n30;
	    const T ds3=-dn1+dn2-n10+n20;

	    
	    const T nsum1=(-dn1-dn1*n20-dn1*n30);
	    //const T nsum1=(n30*n20*(n1+1)-n1*(n20+1)*(n30+1));
	    const T nsum2=Large*(-dn2*n10/Large+dn2*n30/Large);
	    const T nsum3=Large*(-dn3*n10/Large+dn3*n20/Large);
	    

	    Sl[i]+=0.5*dkl*phi*nsum*w[i];
	    dSl[i]+=0.5*dkl*phi*ds1*w[i];

	    //Sl[i]+=0.25*dkl*phi*nsum*w[j2]*dk31/dk32+0.25*dkl*phi*nsum*w3*dk31/dk33;
	    //Sl[i]+=0.5*dkl*phi*nsum*w[j2]+0.5*dkl*phi*nsum*w3;

	    Sl[j2]-=0.5*dkl*phi*nsum*w[j2]*dk31/dk32;
	    dSl[j2]-=0.5*dkl*phi*ds2*w[j2]*dk31/dk32;

	    Sl[j3]-=0.5*dkl*phi*nsum*w3*dk31/dk33;
	    dSl[j3]-=0.5*dkl*phi*ds3*w3*dk31/dk33;

	    //Sl[j2]-=0.25*dkl*phi*nsum*w[j2]*dk31/dk32;
	    //Sl[j3]-=0.25*dkl*phi*nsum*w3*dk31/dk33;//*e012/e0[j3];
	    //S[i]+=0.5*dkl*phi*(n3*n2*(n1+1)-n1*(n2+1)*(n3+1));

	    if(j2==0)
	      int yes(0);
	    if(j3==0)
	      int yes(0);

	  }
      }

    const T preFac=Acell/16./3.141592653*hbarJoule*_maxPhi*hbar*_maxDkl;
    const int klen=_kspace.getlength();
    for(int i=0;i<Rows;i++)
      {
	Sl[i]*=preFac;//*w[i]*hbar);
	dSl[i]*=preFac;
	if(fabs(Sl[i])>fabs(maxS))
	   maxS=fabs(Sl[i]);
      }
    
    long double defect(0);
    int cnt(0);
    T eq(0);
    const int cStart(_kspace.getGlobalIndex(c,0));
    int ind(cStart);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    defect+=((Sl[cnt]/maxS)*(dk3/DK3));
	    //if(fabs(Sl[cnt])*(dk3/DK3)>100)
	    // Sl[cnt]=0.;
	    Tmode& mode=kv.getmode(m);
	    const T tau=mode.gettau();
	    cOt+=mode.calce0(Tl)*(dk3/DK3);
	    //const T e0=mode.calce0(Tl);
	    //eq+=(e0-e[ind])/tau*(dk3/DK3);
	    ind++;
	    cnt++;
	  }
      }

    cOt*=DK3;
    T f(maxS*defect/eq);
    defect*=maxS*DK3;
    //cout<<"Defect1: "<<defect<<endl;
    //cout<<"Temperature: "<<Tl<<endl;
    
    if(correct)
      {
	cnt=0;
	T d(defect);
	defect=0.;
	ind=cStart;
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		const T tau=mode.gettau();
		const T cp=mode.calcde0dT(Tl);
		const T e0(mode.calce0(Tl));
		//S[cnt]+=f*(e0-e[ind])/tau;
		//Sl[cnt]=d*(Sl[cnt]/d-e0/cOt);
		defect+=((Sl[cnt]/maxS)*(dk3/DK3));
		ind++;
		cnt++;
	      }
	  }
	defect*=DK3*maxS;
	//cout<<"Defect2: "<<defect<<endl;
      }

    for(int i=0;i<Rows;i++)
      {
	S[i]=Sl[i];
	dS[i]=dSl[i];
      }

  }

  void IterateToEquilibrium(const T Tl, const int totIts, const T tStep)
  {

    const T hbarJoule=1.054560652926899e-34;
    const T hbar=6.582119e-16;
    const T JouleToeV=6.24150974e18;
    const T Acell=5.378395621705545e-20;
    const T kb=8.617343e-5;  // (eV/K)
    const T cv=1e-9*1e-9;
    T temp=Tl;
    //T Tlat(_kspace.calcLatTemp(0));
    T maxS(0);
    
    TArray S(_kspace.gettotmodes());
    S.zero();
    TArray e0(_kspace.gettotmodes());
    e0.zero();
    TArray e(_kspace.gettotmodes());
    e.zero();
    _kspace.getEquilibriumArray(e,Tl);
    //_kspace.geteCellVals(0,e);
    TArray& w(_kspace.getFreqArray());
    TArray t(_kspace.gettotmodes());
    t.zero();
    TArray gam1(_kspace.gettotmodes());
    gam1.zero();
    TArray gam2(_kspace.gettotmodes());
    gam2.zero();
    TArray gam(_kspace.gettotmodes());
    gam.zero();

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
    const T preFac=1./16.*hbarJoule*_maxDkl*_maxPhi/3.14159265358979323846264;

    for(int it=0;it<totIts;it++)
      {
	S.zero();
	T cOt(0);
	//type 1 collisions
	
	for(int i=0;i<Rows;i++)
	  {
	    for(int pos=t1p2row[i];pos<t1p2row[i+1];pos++)
	      {
		const int j2=t1p2col[pos];
		const int j3=t1p3col[pos];
		const T dkl=t1dkl[pos];
		const T phi=t1phi[pos]*Acell;
		const T w3=w[i]+w[j2];
		const T n1=e[i]/hbar/w[i];
		const T n2=e[j2]/hbar/w[j2];
		const T n3=e[j3]/hbar/w[j3];
		const T n10=1./(exp(hbar*w[i]/kb/300.)-1);
		const T n20=1./(exp(hbar*w[j2]/kb/300.)-1);
		const T n30=1./(exp(hbar*w[j3]/kb/300.)-1);

		S[i]+=dkl*phi*((n1+1)*(n2+1)*n3-n1*n2*(n3+1));
		t[i]+=dkl*phi*(n20-n30);
		gam1[i]+=dkl*phi*n10*n20*(n30+1);
	      }
	  }

	//type 2 collisions
	for(int i=0;i<Rows;i++)
	  {
	    const int k1=floor(i/6);
	    const T dk31=_kspace.getkvol(k1).getdk3();
	    for(int pos=t2p2row[i];pos<t2p2row[i+1];pos++)
	      {
		const int j2=t2p2col[pos];
		const int k2=floor(j2/6);
		const T dk32=_kspace.getkvol(k2).getdk3();
		const int j3=t2p3col[pos];
		const int k3=floor(j3/6);
		const T dk33=_kspace.getkvol(k3).getdk3();
		const T dkl=t2dkl[pos];
		const T phi=t2phi[pos]*Acell;
		const T w3=w[i]-w[j2];
		const T n1=e[i]/w[i]/hbar;
		const T n2=e[j2]/w[j2]/hbar;
		const T n3=e[j3]/w[j3]/hbar;

		/*
		const T n10=e0[i]/w[i]/hbar;
		const T n20=e0[j2]/w[j2]/hbar;
		const T n30=e0[j3]/w[j3]/hbar;
		*/

		const T n10=1./(exp(hbar*w[i]/kb/300.)-1);
		const T n20=1./(exp(hbar*w[j2]/kb/300.)-1);
		const T n30=1./(exp(hbar*w[j3]/kb/300.)-1);

		const T e012=hbar*w3/(exp(hbar*w3/kb/Tl)-1.);
		const T dn1=n1-n10;
		const T dn2=n2-n20;
		const T dn3=(n3-n30);
		

		S[i]+=0.5*dkl*phi*(n3*n2*(n1+1)-n1*(n2+1)*(n3+1));
		t[i]+=0.5*dkl*phi*(n20+n30+1);
		gam2[i]+=0.5*dkl*phi*n10*(n20+1)*(n30+1);
	      }
	  }

	const int klen=_kspace.getlength();
	for(int i=0;i<Rows;i++)
	  {
	    S[i]*=(preFac*w[i]*hbar);
	    t[i]*=preFac;
	    t[i]=1./t[i];
	    const T n10=1./(exp(hbar*w[i]/kb/300.)-1);
	    gam[i]=n10*(n10+1)/preFac/(gam1[i]+gam2[i]);
	    if(fabs(S[i])>fabs(maxS))
	      maxS=fabs(S[i]);
	  }

	const T DK3=_kspace.getDK3();
	T defect(0);
	int cnt(0);
	T maxRatio(0);
	int maxInd(-1);
	int minInd(-1);
	T minRatio(10);
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		const T tau=mode.gettau();
		const T dif=tau-t[cnt];
		const T gamRatio=tau/gam[cnt];
		if(gamRatio>maxRatio)
		  {
		    maxRatio=gamRatio;
		    maxInd=cnt;
		  }
		if(gamRatio<minRatio)
		  {
		    minRatio=gamRatio;
		    minInd=cnt;
		  }
		defect+=(S[cnt]/maxS)*(dk3/DK3);
		cnt++;
	      }
	  }

	cout<<"Max difference: "<<maxRatio<<" at: "<<maxInd<<endl;
	cout<<"Min difference: "<<minRatio<<" at: "<<minInd<<endl;
	
	T eq(0);
	if(it>-1)
	  {
	    cnt=0;
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kv=_kspace.getkvol(k);
		const int modenum=kv.getmodenum();
		T dk3=kv.getdk3();
		for(int m=0;m<modenum;m++)
		  {
		    Tmode& mode=kv.getmode(m);
		    const T tau=mode.gettau();
		    const T cp=mode.calcde0dT(Tl);
		    cOt+=cp/tau*dk3/DK3;
		    eq+=(mode.calce0(temp)-e[cnt])*(dk3/tau/DK3);
		    cnt++;
		  }
	      }

	    defect=defect*maxS*DK3;
	    cOt*=DK3;
	    cout<<"Defect1: "<<defect<<endl;
	    cnt=0;
	    T d(defect);
	    defect=0.;
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kv=_kspace.getkvol(k);
		const int modenum=kv.getmodenum();
		T dk3=kv.getdk3();
		for(int m=0;m<modenum;m++)
		  {
		    Tmode& mode=kv.getmode(m);
		    const T tau=mode.gettau();
		    const T cp=mode.calcde0dT(Tl);
		    S[cnt]=d*(S[cnt]/d-cp/tau/cOt);
		    defect+=(S[cnt]/maxS)*(dk3/DK3);
		    cnt++;
		  }
	      }
	    defect*=DK3*maxS;
	    cout<<"Defect2: "<<defect<<endl;
	  }

	for(int i=0;i<Rows;i++)
	  e[i]+=S[i]*tStep;

	cnt=0;
	T esum(0);
	
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		esum+=e[cnt]*dk3;
		cnt++;
	      }
	  }

	_kspace.calcTemp(temp,esum);
	cout<<"New Temp: "<<temp<<endl;
      }
  }

  void correctDetailedBalance()
  {
    const IntArray& t1p2row=_type1Collisions.getSelfRow();
    const IntArray& t1p3row=_type1Collisions.getOtherRow();
    const IntArray& t1p2col=_type1Collisions.getSelfCol();
    const IntArray& t1p3col=_type1Collisions.getOtherCol();

    const TArray& t1dkl=_type1Collisions.getSelfCoeffs();
    TArray& t1phi=_type1Collisions.getNonConstOtherCoeffs();
    BArray Dups(t1phi.getLength());
    Dups=false;

    const int Rows(_kspace.gettotmodes());
    int duplicated(0);
    int removed(0);

    for(int i=0;i<Rows;i++)
      {
	for(int pos=t1p2row[i];pos<t1p2row[i+1];pos++)
	  {
	    const int j2=t1p2col[pos];
	    const int j3=t1p3col[pos];
	    int p1(-1);
	    int p2(-1);
	    int p3(-1);

	    //type1 for j2
	    bool hasI1(false);
	    if(true)
	      {
		for(int pos2=t1p2row[j2];pos2<t1p2row[j2+1];pos2++)
		  {
		    const int jj2=t1p2col[pos2];
		    if(i==jj2)
		      {
			hasI1=true;
			p1=pos2;
			break;
		      }
		  }
	      }

	    if((hasI1))
	      {
		Dups[pos]=true;
		Dups[p1]=true;
		t1phi[pos]*=0.5;
		duplicated++;
	      }
	  }
      }

    for(int i=0;i<Dups.getLength();i++)
      {
	if(!Dups[i])
	  {
	    //t1phi[i]=0;
	    removed++;
	  }
      }

    cout<<"Duplicated -- "<<duplicated<<endl;
    cout<<"Removed -- "<<removed<<endl;

  }

 private:

  void normalize()
  {
    _type1Collisions.multiplySelf(1./_maxDkl);
    _type1Collisions.multiplyOther(1./_maxPhi);
    _type2Collisions.multiplySelf(1./_maxDkl);
    _type2Collisions.multiplyOther(1./_maxPhi);
  }

  ScatteringKernel(const ScatteringKernel&);
  Kspace<T>& _kspace;
  KSConnectivity<T> _type1Collisions;
  KSConnectivity<T> _type2Collisions;
  T _maxPhi;
  T _maxDkl;

};

#endif
