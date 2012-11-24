#ifndef _SCATTERINGKERNEL_H_
#define _SCATTERINGKERNEL_H_

#include "Array.h"
#include "KSConnectivity.h"
#include <iostream>
#include <fstream>
#include <stdio.h>


template<class T> class Kspace;

template <class T>
class ScatteringKernel
{
  
  typedef Array<T> TArray;
  typedef Array<int> IntArray;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;

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
	
	const T n1=1./(exp(hbar*w[i]/kb/300.)-1);
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    const T err=fabs((n1+1)*(n2+1)*n3-n1*n2*(n3+1))/n1;
	    
	    if(err<tol)
	      {
		_type1Collisions.addCountSelf(i,1);
		_type1Collisions.addCountOther(i,1);
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

	const T n1=1./(exp(hbar*w[i]/kb/300.)-1);
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    const T err=fabs((n1+1)*(n2+1)*n3-n1*n2*(n3+1))/n1;

	      if(err<tol)
		{
		  _type1Collisions.addSelf(i,index2,dkl);
		  _type1Collisions.addOther(i,index3,phi);
		}
	  }
      }

    _type1Collisions.finishAddSelf();
    _type1Collisions.finishAddOther();
    
    cout<<"Type I complete."<<endl;

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
    const int rowLen=_type1Collisions.getSelfSize();
    
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
	
	const T n1=1./(exp(hbar*w[i]/kb/300.)-1);
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    const T err=fabs(n3*n2*(n1+1)-n1*(n2+1)*(n3+1))/n1;

	    if(err<tol)
	      {
		_type2Collisions.addCountSelf(i,1);
		_type2Collisions.addCountOther(i,1);
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

	const T n1=1./(exp(hbar*w[i]/kb/300.)-1);
	//if(i % 2000==0)
	//cout<<"Row: "<<i<<endl;

	fread(&row2,sizeof(int),1,fp_in2);
	fread(&nnz2,sizeof(int),1,fp_in2);
	fread(&row3,sizeof(int),1,fp_in3);
	fread(&nnz3,sizeof(int),1,fp_in3);

	for(int j=0;j<nnz2;j++)
	  {
	    fread(&index2,sizeof(int),1,fp_in2);
	    fread(&dkl,sizeof(T),1,fp_in2);
	    fread(&index3,sizeof(int),1,fp_in3);
	    fread(&phi,sizeof(T),1,fp_in3);

	    const T n2=1./(exp(hbar*w[index2]/kb/300.)-1);
	    const T n3=1./(exp(hbar*w[index3]/kb/300.)-1);
	    const T err=fabs(n3*n2*(n1+1)-n1*(n2+1)*(n3+1))/n1;

	    if(err<tol)
	      {
		_type2Collisions.addSelf(i,index2,dkl);
		_type2Collisions.addOther(i,index3,phi);
	      }
	  }
      }

    _type2Collisions.finishAddSelf();
    _type2Collisions.finishAddOther();
    
    cout<<"Type II complete."<<endl;

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
	    phi=phi/w[i]/w[j2]/w[j3]*pow(hbarJoule,3)/48.;
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
	    phi=phi/w[i]/w[j2]/w[j3]*pow(hbarJoule,3)/48.;
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

  void updateSource(const TArray& e, const TArray& w, TArray& S, const T Tl)
  {

    const T hbarJoule=1.054571726e-34;
    const T hbar=6.582119e-16;
    const T Acell=5.378395621705545e-20;
    T temp=Tl;
    
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
	    const T w3=w[i]+w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w3/hbar;
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
	    const T phi=t2phi[pos];
	    const T w3=w[i]+w[j2];
	    const T n1=e[i]/w[i]/hbar;
	    const T n2=e[j2]/w[j2]/hbar;
	    const T n3=e[j3]/w3/hbar;
	    S[i]+=0.5*dkl*phi*(n3*n2*(n1+1)-n1*(n2+1)*(n3+1));
	  }
      }

    const T preFac=Acell/2./3.141592653/pow(hbarJoule,2)*_maxDkl*_maxPhi;
    const int klen=_kspace.getlength();

    const T DK3=_kspace.getDK3();
    T defect(0);
    int cnt(0);
    T eq(0);
    for(int k=0;k<klen;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    S[cnt]*=preFac*w[cnt]*hbar;
	    defect+=S[cnt]*(dk3/DK3);
	    Tmode& mode=kv.getmode(m);
	    const T tau=mode.gettau();
	    eq+=(mode.calce0(Tl)-e[cnt])*(dk3/tau/DK3);
	    cnt++;
	  }
      }

    eq*=DK3;
    defect*=DK3;

    if(eq>1e14)
      {
	T f(defect/eq);
	cnt=0;
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
		S[cnt]-=f*(mode.calce0(Tl)-e[cnt])/tau;
		defect+=S[cnt]*(dk3/DK3);
		cnt++;
	      }
	  }
	defect*=DK3;
      }
    else
      {
	S.zero();
      }

  }

  void IterateToEquilibrium(const T Tl, const int totIts, const T tStep)
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

    for(int it=0;it<totIts;it++)
      {
	S.zero();

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
		const T w3=w[i]+w[j2];
		const T n1=e[i]/w[i]/hbar;
		const T n2=e[j2]/w[j2]/hbar;
		const T n3=e[j3]/w[j3]/hbar;
		S[i]+=0.5*dkl*phi*(n3*n2*(n1+1)-n1*(n2+1)*(n3+1));
	      }
	  }

	const T preFac=1./2./3.141592653/pow(hbarJoule,2)*_maxDkl*_maxPhi;
	for(int i=0;i<Rows;i++)
	  S[i]*=preFac*w[i]*hbar;

	const int klen=_kspace.getlength();
	const T DK3=_kspace.getDK3();
	T defect(0);
	int cnt(0);
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    for(int m=0;m<modenum;m++)
	      {
		defect+=S[cnt]*(dk3/DK3);
		cnt++;
	      }
	  }

	
	T eq(0);
	if(it>0)
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
		    eq+=(mode.calce0(temp)-e[cnt])*(dk3/tau/DK3);
		    cnt++;
		  }
	      }

	    const T f(defect/eq);
	    cout<<"Factor: "<<f<<endl;
	    cnt=0;
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
		    S[cnt]-=f*(mode.calce0(temp)-e[cnt])/tau;
		    defect+=S[cnt]*(dk3/DK3);
		    cnt++;
		  }
	      }
	    defect*=DK3;

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
	cout<<"Defect: "<<defect<<endl<<endl;
      }
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
