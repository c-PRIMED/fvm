// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KSPACE_H_
#define _KSPACE_H_

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include<string>
#include <vector>
#include <math.h>
#include "Vector.h"
#include "pmode.h"
#include "kvol.h"
#include "SquareTensor.h"
#include "ScatteringKernel.h"

template<class T>
class DensityOfStates;

template<class T>
class Kspace
{

 public:

  typedef Kspace<T> Tkspace;
  typedef vector<Tkspace*> TkspList;
  typedef Array<T> TArray;
  typedef shared_ptr<TArray> TArrPtr;
  typedef Vector<T,3> Tvec;
  typedef Array<Tvec> TvecArray;
  typedef Array<int> IntArray;
  typedef pmode<T> Tmode;
  typedef shared_ptr<Tmode> Tmodeptr;
  typedef vector<Tmodeptr> Modes;
  typedef kvol<T> Tkvol;
  typedef shared_ptr<Tkvol> Kvolptr;
  typedef vector<Kvolptr> Volvec;
  typedef SquareTensor<T,3> T3Tensor;
  typedef typename Tmode::Reflection Reflection;
  typedef typename Tmode::Reflptr Reflptr;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef typename Tmode::Refl_Map Refl_Map;
  typedef pair<TArray*, TArray*> BinAndTrans;
  typedef map<Tkspace*,pair<TArray*,TArray*> > TransmissionMap;
  typedef typename Tkspace::TransmissionMap::iterator TransIt;
  typedef pair<const StorageSite*, const StorageSite*> EntryIndex;
  typedef map<EntryIndex, shared_ptr<ArrayBase> > GhostArrayMap;

 Kspace(T a, T tau, T vgmag, T omega, int ntheta, int nphi) :
  _length(ntheta*nphi),
    _Kmesh(),
    _totvol(0.),
    _freqArray(0)
      { //makes gray, isotropic kspace  
	
	const long double pi=3.141592653589793238462643383279502884197169399;
	long double sumX(0);
	long double sumY(0);
	long double sumZ(0);
	long double theta;
	long double phi;
	long double dtheta=pi/T(ntheta)/2.;
	long double dphi=2.*pi/T(nphi);
	long double dk3;
	const long double Kmax=pi/a;
	const long double Ktot=pi*pow(Kmax,3.)*4./3./pow((2.*pi),3.);
	Tvec vg;
	int count=1;
	for(int t=0;t<ntheta;t++)
	  {
	    theta=dtheta*(T(t)+.5);
	    for(int p=0;p<nphi;p++)
	      {
		phi=dphi*(T(p)+.5);
		vg[0]=vgmag*sin(theta)*sin(phi);
		vg[1]=vgmag*sin(theta)*cos(phi);
		vg[2]=vgmag*cos(theta);
		dk3=(sin(theta)*sin(dtheta/2.)*(dtheta/pi))*Ktot;
		sumX+=vgmag*sin(theta)*sin(phi);
		sumY+=vgmag*sin(theta)*cos(phi);
		sumZ+=vgmag*cos(theta);
		Tmodeptr modeptr=shared_ptr<Tmode>(new Tmode(vg,omega,tau));
		modeptr->setIndex(count);
		count++;
		Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeptr,dk3));
		_Kmesh.push_back(volptr);
		_totvol+=dk3;
	      }
	  }
	cout<<"SumX: "<<sumX<<endl;
	cout<<"SumY: "<<sumY<<endl;
	cout<<"SumZ: "<<sumZ<<endl;
      }

 Kspace():
  _freqArray(0)
    {}

 void setCp(const T cp)
 {//input the total specific heat in eV/m^3/K
   
   for(int k=0;k<_length;k++)
     {
       Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    mode.getcpRef()=cp/_totvol;
	  }
     }
 }

 void setCpNonGray(const T Tl)
 {
   for(int k=0;k<_length;k++)
     {
       Tkvol& kv=getkvol(k);
       const int modenum=kv.getmodenum();
       for(int m=0;m<modenum;m++)
	 {
	   Tmode& mode=kv.getmode(m);
	   const T omega=mode.getomega();
	   const T hbar=6.582119e-16;  // (eV s)
	   const T kb=8.617343e-5;  // (eV/K) 

	   mode.getcpRef()=kb*pow((hbar*omega/kb/Tl),2)*
	     exp(hbar*omega/kb/Tl)/pow((exp(hbar*omega/kb/Tl)-1),2);
	 }
     }
 }

 void makeFreqArray()
 {
   _freqArray.resize(gettotmodes());
   int count=0;
   for(int k=0;k<_length;k++)
     {
       Tkvol& kv=getkvol(k);
       const int modenum=kv.getmodenum();
       for(int m=0;m<modenum;m++)
	 {
	   Tmode& mode=kv.getmode(m);
	   const T omega=mode.getomega();
	   _freqArray[count]=omega;
	   count++;
	 }
     }
 }

 Kspace(const char* filename,const int dimension):
 _freqArray(0)
   {
     ifstream fp_in;
     fp_in.open(filename,ifstream::in);

     int modeNum;
     int kPoints;
     int directions;
     
     cout<<endl<<"Reading BZ file"<<endl;
     fp_in>>modeNum;
     cout<<"Number of Polarizations: "<<modeNum<<endl;
     fp_in>>kPoints;
     cout<<"Number of Wave Vector Magnitude Discretizations: "<<kPoints<<endl;
     fp_in>>directions;
     cout<<"Number of Angular Discretizations: "<<directions<<endl;
     cout<<"Total Number of K-Space Points: "<<modeNum*kPoints*directions<<endl;

     _length=kPoints*directions;
     int count=1;

     for(int k=0;k<_length;k++)
       {
	 Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeNum));
	 Modes& modes=volptr->getModes();

	 for(int m=0;m<modeNum;m++)
	   {
	     T Tdummy(0.);
	     T omega(0.);
	     T tau(0.);
	     Tvec vg;
	     Tvec K; 
	     T weight(0.);
	     Tmodeptr modeptr=shared_ptr<Tmode>(new Tmode());
	     fp_in>>Tdummy;
	     fp_in>>weight;
	     fp_in>>omega;
	     fp_in>>Tdummy;
	     K[0]=Tdummy;

	     if(dimension==2)
	       {
		 fp_in>>Tdummy;
		 K[1]=Tdummy;
		 K[2]=0.;
	       }
	     
	     if(dimension==3)
	       {
		 fp_in>>Tdummy;
		 K[1]=Tdummy;
		 fp_in>>Tdummy;
		 K[2]=Tdummy;
	       }

	     fp_in>>Tdummy;
	     vg[0]=Tdummy;

	     if(dimension==2)
	       {
		 fp_in>>Tdummy;
		 vg[1]=Tdummy;
		 vg[2]=0.;
	       }
	     
	     if(dimension==3)
	       {
		 fp_in>>Tdummy;
		 vg[1]=Tdummy;
		 fp_in>>Tdummy;
		 vg[2]=Tdummy;
	       }

	     fp_in>>tau;

	     modeptr->getVRef()=vg;
	     modeptr->getTauRef()=tau;
	     modeptr->getOmegaRef()=omega;
	     modeptr->setIndex(count);
	     count++;
	     modes.push_back(modeptr);
	     volptr->setkvec(K);
	     volptr->setdk3(weight);
	   }
	 _Kmesh.push_back(volptr);
       }

     fp_in.close();
     calcDK3();
   }
 
 Kspace(const char* filename,const int dimension,const bool normal):
 _freqArray(0)
   {
     ifstream fp_in;
     fp_in.open(filename,ifstream::in);

     int modeNum;
     int kPoints;
     int directions;
     
     cout<<endl<<"Using Shifted Normal Scattering"<<endl;
     cout<<"Reading BZ file"<<endl;
     fp_in>>modeNum;
     cout<<"Number of Polarizations: "<<modeNum<<endl;
     fp_in>>kPoints;
     cout<<"Number of Wave Vector Magnitude Discretizations: "<<kPoints<<endl;
     fp_in>>directions;
     cout<<"Number of Angular Discretizations: "<<directions<<endl;
     cout<<"Total Number of K-Space Points: "<<modeNum*kPoints*directions<<endl;

     _length=kPoints*directions;
     int count=1;

     for(int k=0;k<_length;k++)
       {
	 Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeNum));
	 Modes& modes=volptr->getModes();

	 for(int m=0;m<modeNum;m++)
	   {
	     T Tdummy(0.);
	     T omega(0.);
	     T tau(0.);
	     T tauN(0.);
	     Tvec vg;
	     Tvec K; 
	     T weight(0.);
	     Tmodeptr modeptr=shared_ptr<Tmode>(new Tmode());
	     fp_in>>Tdummy;
	     fp_in>>weight;
	     fp_in>>omega;
	     fp_in>>Tdummy;
	     K[0]=Tdummy;

	     if(dimension==2)
	       {
		 fp_in>>Tdummy;
		 K[1]=Tdummy;
		 K[2]=0.;
	       }
	     
	     if(dimension==3)
	       {
		 fp_in>>Tdummy;
		 K[1]=Tdummy;
		 fp_in>>Tdummy;
		 K[2]=Tdummy;
	       }

	     fp_in>>Tdummy;
	     vg[0]=Tdummy;

	     if(dimension==2)
	       {
		 fp_in>>Tdummy;
		 vg[1]=Tdummy;
		 vg[2]=0.;
	       }
	     
	     if(dimension==3)
	       {
		 fp_in>>Tdummy;
		 vg[1]=Tdummy;
		 fp_in>>Tdummy;
		 vg[2]=Tdummy;
	       }

	     fp_in>>tau;
	     fp_in>>tauN;

	     modeptr->getVRef()=vg;
	     modeptr->getTauRef()=tau;
	     modeptr->getTauNRef()=tauN;
	     modeptr->getOmegaRef()=omega;
	     modeptr->setIndex(count);
	     count++;
	     modes.push_back(modeptr);
	     volptr->setkvec(K);
	     volptr->setdk3(weight);
	   }
	 _Kmesh.push_back(volptr);
       }

     fp_in.close();
     calcDK3();
   }
  
  //void setvol(int n,Tkvol k) {*_Kmesh[n]=k;}
  Tkvol& getkvol(int n) const {return *_Kmesh[n];}
  int getlength() const {return _length;}
  T gethbar() {return 6.582119e-16;}
  int gettotmodes()
  {
    return (_Kmesh[0]->getmodenum())*_length;
  }
  T calcDK3()
  {
    T r(0.0);
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	r+=kv.getdk3();
      }
    _totvol=r;
    return r;
  }
  T getDK3() const {return _totvol;}
  T calcTauTot()
  {   // returns sum(dk3/tau)
    T tauTot=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    tauTot+=dk3/mode.gettau();
	  }
      }
    return tauTot;
  }
  void NewtonSolve(T& guess, const T e_sum)
  {
    T e0_tau;
    T de0_taudT;
    T deltaT=1.;
    T newguess;
    int iters=0;
    while((deltaT>1e-6)&&(iters<10))
      {
	gete0_tau(guess,e0_tau,de0_taudT);
	deltaT=(e0_tau-e_sum)/de0_taudT;
	newguess=guess-deltaT;
	deltaT=fabs(deltaT/guess);
	guess=newguess;
	iters++;
      }	
  }

  void calcTemp(T& guess, const T e_sum, const Tvec An)
  {
    T e0;
    T de0dT;
    T deltaT=1.;
    T newguess;
    int iters=0;
    while((deltaT>1e-6)&&(iters<10))
      {
	gete0(guess, e0, de0dT,An);
	deltaT=(e0-e_sum)/de0dT;
	newguess=guess-deltaT;
	deltaT=fabs(deltaT/guess);
	guess=newguess;
	iters++;
      }	
  }

  T calcModeTemp(T guess, const T e_sum, const T m)
  {
    T e0;
    T de0dT;
    T deltaT=1.;
    T newguess;
    int iters=0;
    while((deltaT>1e-6)&&(iters<10))
      {
	gete0(guess, e0, de0dT,m);
	deltaT=(e0-e_sum)/de0dT;
	newguess=guess-deltaT;
	deltaT=fabs(deltaT/guess);
	guess=newguess;
	iters++;
      }
    return guess;
  }

  void gete0_tau(T& Tguess, T& e0tau, T& de0taudT)
  {
    e0tau=0.;
    de0taudT=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    e0tau+=mode.calce0tau(Tguess)*dk3;
	    de0taudT+=mode.calcde0taudT(Tguess)*dk3;
	  }
      }
  }

  void gete0(const T Tguess, T& e0, T& de0dT, const Tvec An)
  {
    e0=0.;
    de0dT=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    Tvec vg=mode.getv();
	    T vdota=vg[0]*An[0]+vg[1]*An[1]+vg[2]*An[2];
	    if(vdota>0)
	      {
		e0+=mode.calce0(Tguess)*dk3;
		de0dT+=mode.calcde0dT(Tguess)*dk3;
	      }
	  }
      }
  }

  void gete0(const T Tguess, T& e0, T& de0dT, const T m)
  {
    e0=0.;
    de0dT=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	T dk3=kv.getdk3();
	Tmode& mode=kv.getmode(m);
	e0+=mode.calce0(Tguess)*dk3;
	de0dT+=mode.calcde0dT(Tguess)*dk3;
      }
  }

  T getde0taudT(T Tl)
  {
    T de0taudT=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    de0taudT+=mode.calcde0taudT(Tl)*dk3;
	  }
      }
    return de0taudT;
  }

  /*
  T getde0taudTgray()
  {
    T de0taudT=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    de0taudT+=mode.calcde0taudTgray()*dk3;
	  }
      }
    return de0taudT;
  }
  */

  T calcSpecificHeat(T Tl)
  {
    T r(0.0);
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    r+=mode.calcde0dT(Tl)*dk3;
	  }
      }
    return r;
  }

  T calcSpecificHeat(T Tl,const int m)
  {
    T r(0.0);
    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	T dk3=kv.getdk3();
	Tmode& mode=kv.getmode(m);
	r+=mode.calcde0dT(Tl)*dk3;
      }
    return r;
  }

  T findKnStats(const T length)
  {
    T AveKn(0.0);
    T maxKn(0.0);
    T minKn;

    Tvec vg1=getkvol(0).getmode(0).getv();
    T vmag1=sqrt(pow(vg1[0],2)+pow(vg1[1],2)+pow(vg1[2],2));
    T tau1=getkvol(0).getmode(0).gettau();
    T npol(getkvol(0).getmodenum());
    minKn=vmag1*tau1;
    if(minKn==0.)
      minKn=1.;

    for(int k=0;k<_length;k++)
      {
	Tkvol& kv=getkvol(k);
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	for(int m=0;m<modenum;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    Tvec vg=mode.getv();
	    T vmag=sqrt(pow(vg[0],2)+pow(vg[1],2)+pow(vg[2],2));
	    T tau=mode.gettau();
	    AveKn+=vmag*tau*dk3/_totvol;
	    if(vmag*tau>maxKn)
	      maxKn=vmag*tau;
	    if(vmag*tau<minKn && vmag*tau>0.)
	      minKn=vmag*tau;

	  }
      }
    AveKn/=length;
    AveKn/=npol;
    maxKn/=length;
    minKn/=length;

    cout<<"Average Kn: "<<AveKn<<endl;
    cout<<"Maximum Kn: "<<maxKn<<endl;
    cout<<"Minimum Kn: "<<minKn<<endl;

    return AveKn;
    
  }

  void findSpecs(T dk3, T vo, int m, Tvec so, Refl_pair& refls)
  {
    Tmode& mode1=getkvol(0).getmode(m);
    Tmode& mode2=getkvol(1).getmode(m);
    int m1=0;
    int m2=1;
    Tvec vg1=mode1.getv();
    Tvec vg2=mode2.getv();
    Tvec sn1=vg1;///sqrt(pow(vg1[0],2)+pow(vg1[1],2)+pow(vg1[2],2));
    Tvec sn2=vg2;///sqrt(pow(vg2[0],2)+pow(vg2[1],2)+pow(vg2[2],2));
    T sn1dotso=fabs(sn1[0]*so[0]+sn1[1]*so[1]+sn1[2]*so[2]-pow(vo,2));
    T sn2dotso=fabs(sn2[0]*so[0]+sn2[1]*so[1]+sn2[2]*so[2]-pow(vo,2));
    
    if (sn2dotso<sn1dotso)
      {
	T temp=sn1dotso;
	sn1dotso=sn2dotso;
	sn2dotso=temp;
	m1=1;
	m2=0;
      }

    for(int k=2;k<_length;k++)
      {
	Tmode& temp_mode=getkvol(k).getmode(m);
	Tvec vgt=temp_mode.getv();
	const Tvec sn=vgt;///sqrt(pow(vgt[0],2)+pow(vgt[1],2)+pow(vgt[2],2));
	const T sndotso=fabs(sn[0]*so[0]+sn[1]*so[1]+sn[2]*so[2]-pow(vo,2));
	
	if(sndotso<sn1dotso)
	  {
	    sn2dotso=sn1dotso;
	    sn1dotso=sndotso;
	    m2=m1;
	    m1=k;
	  }
	else if (sndotso<sn2dotso)
	  {
	    sn2dotso=sndotso;
	    m2=k;
	  }
      }
    
    T w1=sn1dotso/(sn1dotso+sn2dotso);
    T w2=sn2dotso/(sn1dotso+sn2dotso);

    //cout<<"Closest: "<<m1<<","<<sn1dotso<<endl;
    
    if(sn1dotso>.99)
      {
	w1=1;
	w2=0;
      }

    T dk31=getkvol(m1).getdk3();
    T dk32=getkvol(m2).getdk3();

    vg1=getkvol(m1).getmode(m).getv();
    vg2=getkvol(m2).getmode(m).getv();

    T v1mag=sqrt(pow(vg1[0],2)+pow(vg1[1],2)+pow(vg1[2],2));
    T v2mag=sqrt(pow(vg2[0],2)+pow(vg2[1],2)+pow(vg2[2],2));

    refls.first.first=w1*vo*dk3/v1mag/dk31;
    refls.second.first=w2*vo*dk3/v2mag/dk32;
    refls.first.second=m1; //refls.first-- to whom the mode dumps energy
    refls.second.second=-1;  //refls.second-- from whom the mode receices energy 
  }
  
  void CopyKspace(Tkspace& copyFromKspace)
  {
    _length=copyFromKspace.getlength();
    _totvol=copyFromKspace.getDK3();
    _Kmesh.clear();
    for(int i=0;i<_length;i++)
      {
	Kvolptr newKvol=Kvolptr(new Tkvol());
	newKvol->copyKvol(copyFromKspace.getkvol(i));
	_Kmesh.push_back(newKvol);
      }

    setDOS(*copyFromKspace.getDOSptr());

  }

  T FindBallisticHeatRate(const Tvec An,const T T1,const T T2)
  {
    T q=0.;
    for(int k=0;k<_length;k++)
      {
	Tkvol& kvol=getkvol(k);
	const T dk3=kvol.getdk3();
	const int modes=kvol.getmodenum();
	for(int m=0;m<modes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Tvec vg=mode.getv();
	    T flux=vg[0]*An[0]+vg[1]*An[1]+vg[2]*An[2];
	    T e0;
	    if(flux>0.)
	      e0=mode.calce0(T2);
	    else
	      e0=mode.calce0(T1);
	    q+=flux*e0*dk3;
	  }
      }
    return q;
  }

  ArrayBase* getVelocities()
  {
    const int allModes=gettotmodes();
    TvecArray* Velocities=new TvecArray(allModes);
 
    for(int k=0;k<_length;k++)
      {
	Tkvol& kvol=getkvol(k);
	const int modes=kvol.getmodenum();
	for(int m=0;m<modes;m++)
	  {
	    const int count=kvol.getmode(m).getIndex()-1;
	    (*Velocities)[count]=kvol.getmode(m).getv();
	  }
      }
    return Velocities;
  }

  ArrayBase* getReflectionArray(const Mesh& mesh, const int FgId)
    {
      const int allModes=gettotmodes();
      IntArray* reflInd=new IntArray(allModes);

      for(int k=0;k<_length;k++)
	{
	  Tkvol& kvol=getkvol(k);
	  const int modes=kvol.getmodenum();
	  for(int m=0;m<modes;m++)
	    {
	      Tmode& mode=kvol.getmode(m);
	      const int count=mode.getIndex()-1;
	      Refl_pair& refls=mode.getReflpair(FgId);
	      if(refls.second.second!=-1)  //v dot A < 0
		{
		  const int kk=refls.second.second;
		  Tmode& FromMode=getkvol(kk).getmode(m);
		  const int indx=FromMode.getIndex();
		  (*reflInd)[count]=indx;
		}
	      else if(refls.first.second!=-1)//v dot A > 0
		{
		  const int kk=refls.first.second;
		  Tmode& ToMode=getkvol(kk).getmode(m);
		  const int indx=ToMode.getIndex();
		  (*reflInd)[count]=indx;
		}
	      else
		throw CException("Not a reflecting wall!");
	    }
	}
      return reflInd;
    }

  ArrayBase* getHollandConductivity(const T Tl)
  {//returns the thermal conductivity tensor in row major order
    T3Tensor KTensor;
    T3Tensor Dummy;
    KTensor.zero();
    Dummy.zero();
    
    for(int k=0;k<_length;k++)
      {
	Tkvol& kvol=getkvol(k);
	const T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Tvec vg=mode.getv();
	    T tau=mode.gettau();
	    T de0dT=mode.calcde0dT(Tl);
	    outerProduct(vg, vg, Dummy);
	    KTensor+=Dummy*tau*de0dT*dk3;
	    Dummy.zero();
	  }
      }

    TArray* Kptr=new TArray(9);
    int count=0;
    for(int j=0;j<3;j++)
      {
	for(int i=0;i<3;i++)
	  {
	    (*Kptr)[count]=KTensor(i,j);
	    count++;
	  }
      }

    return Kptr;
  }

  void outerProduct(const Tvec& v1, const Tvec& v2, T3Tensor& out)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	out(i,j)=v1[i]*v2[j];
  }

  void setTransmission(Tkspace& toKspace, ArrayBase* freqBins, ArrayBase* transArray)
  {
    BinAndTrans* BTPtr=new BinAndTrans;
    BTPtr->first=dynamic_cast<TArray*>(freqBins);
    BTPtr->second=dynamic_cast<TArray*>(transArray);
    _transMap[&toKspace]=*BTPtr;
  }
  
  T findTransmission(Tkspace& toKspace, const T freq)
  {
    TArray& freqBins=*_transMap[&toKspace].first;
    TArray& trans=*_transMap[&toKspace].second;

    int minInd=0;
    int maxInd=freqBins.getLength()-1;

    if(freq<freqBins[minInd] || freq>freqBins[maxInd])
      throw CException("Frequency not in the given range!");

    while(1)
      {
	int mid=floor((minInd+maxInd)/2.);
	T midFreq=freqBins[mid];

	if(freq>midFreq)
	  {
	    minInd=mid;
	    if(freq<=freqBins[minInd+1])
	      return trans[minInd];
	  }
	else
	  {
	    maxInd=mid;
	    if(freq>freqBins[maxInd-1])
	      return trans[maxInd-1];
	  }
      }

  }
  
  TArray& getTransArray(Tkspace& toKspace)
    {return *_transMap[&toKspace].second;}
  
  T calcBallisticInterface(Tkspace& kspace1, const Tvec& An, const T T0, const T T1)
  {
    T heatRate0(0.);
    T heatRate1(0.);
    const int k1len=kspace1.getlength();
    const T DK0=getDK3();
    
    for(int k=0;k<_length;k++)
      {
	Tkvol& kvol=getkvol(k);
	T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Tvec vg=mode.getv();
	    const T VdotA=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
	    
	    if(VdotA>0.)
	      {
		const T e0=mode.calce0(T0);
		const T t01=findTransmission(kspace1,mode.getomega());
		heatRate0+=t01*VdotA*dk3*e0/DK0;
	      }
	  }
      }
    
    for(int k=0;k<k1len;k++)
      {
	Tkvol& kvol=kspace1.getkvol(k);
	T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Tvec vg=mode.getv();
	    const T VdotA=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
	    
	    if(VdotA<0.)
	      {
		const T e0=mode.calce0(T1);
		const T t10=kspace1.findTransmission(*this,mode.getomega());
		heatRate1+=t10*VdotA*dk3*e0/DK0;
	      }
	  }
      }

    heatRate1/=heatRate0;
    heatRate1=1+heatRate1;
    heatRate1*=heatRate0*DK0;
    
    return heatRate1;
  }

  T calcDiffuseE(Tkspace& kspace1, const Tvec& An, const T T0, const T T1)
  {
    T heatRate0(0.);
    T heatRate1(0.);
    const int k1len=kspace1.getlength();
    const T DK0=getDK3();
    
    for(int k=0;k<_length;k++)
      {
	Tkvol& kvol=getkvol(k);
	T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Tvec vg=mode.getv();
	    const T VdotA=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
	    
	    if(VdotA>0.)
	      {
		const T e0=mode.calce0(T0);
		const T t01=findTransmission(kspace1,mode.getomega());
		heatRate0+=t01*VdotA*dk3*e0/DK0;
	      }
	  }
      }
    
    T sumVdotA(0.);

    for(int k=0;k<k1len;k++)
      {
	Tkvol& kvol=kspace1.getkvol(k);
	T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Tvec vg=mode.getv();
	    const T VdotA=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
	    
	    if(VdotA<0.)
	      {
		const T e0=mode.calce0(T1);
		const T t10=kspace1.findTransmission(*this,mode.getomega());
		heatRate1-=(1.-t10)*VdotA*dk3*e0/DK0;
	      }
	    else
	      sumVdotA+=VdotA*dk3/DK0;
	  }
      }

    heatRate1=heatRate0/sumVdotA+heatRate1/sumVdotA;
    
    return heatRate1;
  }

  void giveTransmissions()
  {
    TransmissionMap& coarseTrans=_coarseKspace->getTransMap();

    for(TransIt it=_transMap.begin();it!=_transMap.end();it++)
      {
	Tkspace* fineToKspace=it->first;
	Tkspace* coarseToKspace=fineToKspace->getCoarseKspace();
	coarseTrans[coarseToKspace]=(it->second);
      }
  }

  void setDOS(DensityOfStates<T>& DOS) {_DOS=&DOS;}
  DensityOfStates<T>* getDOSptr() {return _DOS;}
  void setCoarseKspace(Tkspace* cK) {_coarseKspace=cK;}
  Tkspace* getCoarseKspace() {return _coarseKspace;}
  TransmissionMap& getTransMap() {return _transMap;}
  T& gete(const int cell, const int count) {return (*_e)[cell*gettotmodes()+count];}
  T& gete0(const int cell, const int count) {return (*_e0)[cell*gettotmodes()+count];}
  T& getInj(const int cell, const int count) {return (*_injected)[cell*gettotmodes()+count];}
  T& getRes(const int cell, const int count) {return (*_residual)[cell*gettotmodes()+count];}
  T& getFas(const int cell, const int count) {return (*_FASCorrection)[cell*gettotmodes()+count];}
  int getGlobalIndex(const int cell, const int count) {return cell*gettotmodes()+count;}
  void seteArray(TArrPtr ePtr) {_e=ePtr;}
  void sete0Array(TArrPtr e0Ptr) {_e0=e0Ptr;}
  void setInjArray(TArrPtr InjPtr) {_injected=InjPtr;}
  void setResArray(TArrPtr ResPtr) {_residual=ResPtr;}
  void setFASArray(TArrPtr FASPtr) {_FASCorrection=FASPtr;}
  void setTauArray(TArrPtr TauPtr) {_Tau=TauPtr;}
  const TArray& getFreqArray() {return _freqArray;}
  TArray& geteArray() {return *_e;}
  TArray& gete0Array() {return *_e0;}
  TArray& getInjArray() {return *_injected;}
  TArray& getResArray() {return *_residual;}
  TArray& getFASArray() {return *_FASCorrection;}
  void geteCellVals(const int c, TArray& o)
  {
    int start=getGlobalIndex(c,0);
    int end=start+gettotmodes();
    for(int i=start; i<end; i++)
      o[i-start]=(*_e)[i];
  }

  void seteCellVals(const int c, const TArray& o)
  {
    int start=getGlobalIndex(c,0);
    int end=start+gettotmodes();
    for(int i=start; i<end; i++)
      (*_e)[i]=o[i-start];
  }

  void setResidCell(const int c, const TArray& o)
  {
    int start=getGlobalIndex(c,0);
    int end=start+gettotmodes();
    for(int i=start; i<end; i++)
      (*_residual)[i]=o[i-start];
  }

  void addFAS(const int c, TArray& Bvec)
  {
    int start=getGlobalIndex(c,0);
    int end=start+gettotmodes();
    for(int i=start; i<end; i++)
      Bvec[i-start]-=(*_FASCorrection)[i];
  }

  void addFASint(const int c, TArray& Bvec)
  {
    int start=getGlobalIndex(c,0);
    int end=start+gettotmodes();
    for(int i=start; i<end; i++)
      Bvec[i-start]+=(*_FASCorrection)[i];
  }

  void makeFAS() {(*_FASCorrection)-=(*_residual);}

  void syncLocal(const StorageSite& site)
  {
    //package values to be sent/received
    syncScatter(site);
    createSyncGather(site);

#ifdef FVM_PARALLEL
    //SENDING
    MPI::Request   request_send[get_request_size(site)];
    MPI::Request   request_recv[get_request_size(site)];
    int indxSend = 0;
    int indxRecv = 0;

    const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
    foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
      {
	const StorageSite&  oSite = *mpos.first;
	//checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
	EntryIndex e(&site,&oSite);
	ArrayBase& sendArray = *_ghostArrays[e];
	int to_where  = oSite.getGatherProcID();
	if ( to_where != -1 )
	  {
	    int mpi_tag = oSite.getTag();
	    request_send[indxSend++] =  
	      MPI::COMM_WORLD.Isend( sendArray.getData(), sendArray.getDataSize(), MPI::BYTE, to_where, mpi_tag );
	  }
      }

    //RECIEVING
    //getting values from other meshes to fill g
    const StorageSite::GatherMap& gatherMap = site.getGatherMap();
    foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
      {
	const StorageSite&  oSite = *mpos.first;
	//checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
	EntryIndex e(&oSite,&site);
	ArrayBase& recvArray = *_ghostArrays[e];
	int from_where       = oSite.getGatherProcID();
	if ( from_where != -1 )
	  {
	    int mpi_tag = oSite.getTag();
	    request_recv[indxRecv++] =  
	      MPI::COMM_WORLD.Irecv( recvArray.getData(), recvArray.getDataSize(), MPI::BYTE, from_where, mpi_tag );
	  }
      }
    

    int count  = get_request_size(site);
    MPI::Request::Waitall( count, request_recv );
    MPI::Request::Waitall( count, request_send );
#endif

    syncGather(site);

  }

  void syncScatter(const StorageSite& site)
  {
    const StorageSite::ScatterMap& scatterMap = site.getScatterMap();

    const int totModes=gettotmodes();

    foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
      {
	const StorageSite& oSite = *mpos.first;
	EntryIndex e(&site, &oSite);	
	const Array<int>& fromIndices = *(mpos.second);

	if (_ghostArrays.find(e) == _ghostArrays.end())
	  _ghostArrays[e] = _e->newSizedClone( fromIndices.getLength()*totModes);

	const int cellCount = fromIndices.getLength();	 
	TArray& ghostArray = dynamic_cast<TArray&>(*_ghostArrays[e]);
	
	int ghostIndex(0);
	for(int index=0;index<cellCount;index++)
	  {
	    const int c=fromIndices[index];
	    for(int cKindex=getGlobalIndex(c,0);
		cKindex<getGlobalIndex(c,0)+totModes;
		cKindex++)
	      {
		ghostArray[ghostIndex]=(*_e)[cKindex];
		ghostIndex++;
	      }
	  }

      }
  }

  void createSyncGather(const StorageSite& site)
  {
    const StorageSite::GatherMap& gatherMap = site.getGatherMap();

    foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
      {
	const StorageSite& oSite = *mpos.first;
	EntryIndex e(&oSite, &site);
	const Array<int>& toIndices = *(mpos.second);
	if (_ghostArrays.find(e) == _ghostArrays.end())
	  _ghostArrays[e] = _e->newSizedClone(toIndices.getLength()*gettotmodes());
      }

  }
  
  void syncGather(const StorageSite& site)
  {
    const StorageSite::GatherMap& gatherMap = site.getGatherMap();

    const int totModes=gettotmodes();

    foreach(const StorageSite::GatherMap::value_type& mpos, gatherMap)
      {
	const StorageSite& oSite = *mpos.first;
	const Array<int>& toIndices = *(mpos.second);
	EntryIndex e(&oSite, &site);

	if (_ghostArrays.find(e) == _ghostArrays.end())
	  {
	    ostringstream err;
	    err << "Kspace::syncScatter: ghost array not found for"
	      << &oSite << endl;
	    throw CException(err.str());
	  }

	const TArray& ghostArray=dynamic_cast<const TArray&>(*_ghostArrays[e]);
	const int cellCount=toIndices.getLength();

	int ghostIndex(0);
	for(int index=0;index<cellCount;index++)
	  {
	    const int c=toIndices[index];
	    for(int cKindex=getGlobalIndex(c,0);
		cKindex<getGlobalIndex(c,0)+totModes;
		cKindex++)
	      {
		(*_e)[cKindex]=ghostArray[ghostIndex];
		ghostIndex++;
	      }
	  }

      }

  }

  int get_request_size(const StorageSite& site)
  {
    int indx(0);
    const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
    foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap)
      {
	const StorageSite&  oSite = *mpos.first;
	//checking if storage site is only site or ghost site, we only communicate ghost site ( oSite.getCount() == -1 ) 
	if ( oSite.getGatherProcID() != -1 )
	  indx++;
      }
    return indx++;
  }
  
 private:

  Kspace(const Kspace&);
  //num volumes
  int _length;
  Volvec _Kmesh;
  T _totvol;    //total Kspace volume
  TransmissionMap _transMap;
  DensityOfStates<T>* _DOS;
  Tkspace* _coarseKspace;
  TArrPtr _e;
  TArrPtr _e0;
  TArrPtr _injected;
  TArrPtr _residual;
  TArrPtr _FASCorrection;
  TArrPtr _Tau;
  GhostArrayMap _ghostArrays;
  TArray _freqArray;
  
};


#endif
