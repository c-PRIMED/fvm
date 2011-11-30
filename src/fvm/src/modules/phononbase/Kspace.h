#ifndef _KSPACE_H_
#define _KSPACE_H_

#include <iostream>
#include <fstream>
#include<string>
#include <vector>
#include <math.h>
#include "Vector.h"
#include "pmode.h"
#include "kvol.h"

template<class T>
class Kspace
{

 public:

  typedef Kspace<T> Tkspace;
  typedef Vector<T,3> Tvec;
  typedef pmode<T> Tmode;
  typedef shared_ptr<Tmode> Tmodeptr;
  typedef vector<Tmodeptr> Modes;
  typedef kvol<T> Tkvol;
  typedef shared_ptr<Tkvol> Kvolptr;
  typedef vector<Kvolptr> Volvec;
  typedef typename Tmode::Reflection Reflection;
  typedef typename Tmode::Reflptr Reflptr;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef typename Tmode::Refl_Map Refl_Map;

 Kspace(T a, T tau, T vgmag, T omega, int ntheta, int nphi) :
  _length(ntheta*nphi),
    _Kmesh(),
    _totvol(0.)
      { //makes gray, isotropic kspace  
	
	const double pi=3.141592653;
	T theta;
	T phi;
	T dtheta=pi/ntheta;
	T dphi=2.*pi/nphi;
	T dk3;
	Tvec vg;
	for(int t=0;t<ntheta;t++)
	  {
	    theta=dtheta*(t+.5);
	    for(int p=0;p<nphi;p++)
	      {
		phi=dphi*(p+.5);
		vg[0]=vgmag*sin(theta)*sin(phi);
		vg[1]=vgmag*sin(theta)*cos(phi);
		vg[2]=vgmag*cos(theta);
		dk3=2.*sin(theta)*sin(dtheta/2.)*dphi;
		Tmodeptr modeptr=shared_ptr<Tmode>(new Tmode(vg,omega,tau));
		Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeptr,dk3));
		_Kmesh.push_back(volptr);
		_totvol+=dk3;
	      }
	  }
      }

 Kspace()
      {}

 Kspace(const char* filename,const int dimension)
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

     for(int k=0;k<_length;k++)
       {
	 Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeNum));
	 Modes& modes=volptr->getModes();

	 for(int m=0;m<modeNum;m++)
	   {
	     T Tdummy;
	     T omega;
	     T tau;
	     Tvec vg;
	     Tvec K; 
	     T weight;
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
	     modes.push_back(modeptr);
	     volptr->setkvec(K);
	     volptr->setdk3(weight);
	   }
	 _Kmesh.push_back(volptr);
       }

     fp_in.close();
     calcDK3();
   }
 
 Kspace(const char* filename,const int dimension,const bool normal)
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

     for(int k=0;k<_length;k++)
       {
	 Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeNum));
	 Modes& modes=volptr->getModes();

	 for(int m=0;m<modeNum;m++)
	   {
	     T Tdummy;
	     T omega;
	     T tau;
	     T tauN;
	     Tvec vg;
	     Tvec K; 
	     T weight;
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
    while((deltaT>1e-6)&&(iters<100))
      {
	gete0_tau(guess,e0_tau,de0_taudT);
	deltaT=(e0_tau-e_sum)/de0_taudT;
	newguess=guess-deltaT;
	deltaT=fabs(deltaT/guess);
	guess=newguess;
	iters++;
      }	
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
	const int modenum=kv.getmodenum();
	T dk3=kv.getdk3();
	Tmode& mode=kv.getmode(m);
	r+=mode.calcde0dT(Tl)*dk3;
      }
    return r;
  }

  void findKnStats(const T length)
  {
    T AveKn(0.0);
    T maxKn(0.0);
    T minKn;

    Tvec vg1=getkvol(0).getmode(0).getv();
    T vmag1=sqrt(pow(vg1[0],2)+pow(vg1[1],2)+pow(vg1[2],2));
    T tau1=getkvol(0).getmode(0).gettau();
    minKn=vmag1*tau1;

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
	    AveKn+=vmag*tau*dk3;
	    if(vmag*tau>maxKn)
	      maxKn=vmag*tau;
	    if(vmag*tau<minKn)
	      minKn=vmag*tau;

	  }
      }
    AveKn/=_totvol;
    AveKn/=length;
    maxKn/=length;
    minKn/=length;

    cout<<"Average Kn: "<<AveKn<<endl;
    cout<<"Maximum Kn: "<<maxKn<<endl;
    cout<<"Minimum Kn: "<<minKn<<endl;
    
  }

  void findSpecs(T dk3, T vo, int m, Tvec so, Refl_pair& refls)
  {
    Tmode& mode1=getkvol(0).getmode(m);
    Tmode& mode2=getkvol(1).getmode(m);
    int m1=0;
    int m2=1;
    Tvec vg1=mode1.getv();
    Tvec vg2=mode2.getv();
    Tvec sn1=vg1/sqrt(pow(vg1[0],2)+pow(vg1[1],2)+pow(vg1[2],2));
    Tvec sn2=vg2/sqrt(pow(vg2[0],2)+pow(vg2[1],2)+pow(vg2[2],2));
    T sn1dotso=sn1[0]*so[0]+sn1[1]*so[1]+sn1[2]*so[2];
    T sn2dotso=sn2[0]*so[0]+sn2[1]*so[1]+sn2[2]*so[2];
    
    if (sn2dotso>sn1dotso)
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
	const Tvec sn=vgt/sqrt(pow(vgt[0],2)+pow(vgt[1],2)+pow(vgt[2],2));
	const T sndotso=sn[0]*so[0]+sn[1]*so[1]+sn[2]*so[2];
	
	if(sndotso>sn1dotso)
	  {
	    sn2dotso=sn1dotso;
	    sn1dotso=sndotso;
	    m2=m1;
	    m1=k;
	  }
	else if (sndotso>sn2dotso)
	  {
	    sn2dotso=sndotso;
	    m2=k;
	  }
      }
    
    T w1=sn1dotso/(sn1dotso+sn2dotso);
    T w2=sn2dotso/(sn1dotso+sn2dotso);

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
    refls.first.second=m1;
    refls.second.second=m2;  
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
  }
  
 private:

  Kspace(const Kspace&);
  //num volumes
  int _length;
  Volvec _Kmesh;
  T _totvol;    //total Kspace volume
  
};


#endif
