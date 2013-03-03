// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PMODE_H_
#define _PMODE_H_

#include "Vector.h"
#include "Field.h"
#include <map>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

template<class T>
class pmode
{
  
 public:
  
  typedef pmode<T> Tmode;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Vector<T_Scalar,3> Tvec;
  typedef shared_ptr<pmode<T> > Mode_ptr;
  typedef pair<T,int> Reflection;  //<weight,kvol index>
  typedef shared_ptr<Reflection> Reflptr;
  typedef pair<Reflection,Reflection> Refl_pair;  //two nearest vectors
  typedef map<int,Refl_pair> Refl_Map;  //list mapped to each fg.id
  
 pmode(Tvec vg,T omega,T tau):
  _vg(vg),
    _tau(tau),
    _omega(omega),
    _efield("edoubleprime"),
    _eShifted("eShifted"),
    _e0field("e0"),
    _injected("injected"),
    _residual("residual"),
    _FASCorrection("FASCorrection"),
    _reflections(),
    _index(),
    _Tref(299.)
      {}
    
 pmode():
  _efield("edoubleprime"),
    _eShifted("eShifted"),
    _e0field("e0"),
    _injected("injected"),
    _residual("residual"),
    _FASCorrection("FASCorrection"),
    _reflections(),
    _index(),
    _Tref(299.)
      {};
      
  Tvec getv(){return _vg;}
  T getcp() {return _cp;}
  T gettau() {return _tau;}
  T gettauN() {return _tauN;}
  T getomega() {return _omega;}
  Tvec& getVRef(){return _vg;}
  T& getTauRef() {return _tau;}
  T& getTauNRef() {return _tauN;}
  T& getOmegaRef() {return _omega;}
  T& getcpRef() {return _cp;}
  Refl_Map& getreflmap() {return _reflections;}
  Refl_Map getreflmapValue() {return _reflections;}
  Refl_pair& getReflpair(int i) {return _reflections[i];}
  void setIndex(int index) {_index=index;}
  int getIndex() {return _index;}
  Field& getfield() {return _efield;}
  Field& geteShifted() {return _eShifted;}
  Field& gete0field() {return _e0field;}
  Field& getInjected() {return _injected;}
  Field& getresid() {return _residual;}
  Field& getFASfield() {return _FASCorrection;}
  void setTref(const T Tref) {_Tref=Tref;}

  //========================//
  //  Non Constant Cp Methods     //
  //========================//

      
  /*
    T calce0(T Tl)
    {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K)
    T e0kp;
	
    e0kp=hbar*_omega/(exp(hbar*_omega/kb/Tl)-1);
    //e0kp=1./(exp(hbar*_omega/kb/Tl)-1);
    //T eref=1./(exp(hbar*_omega/kb/300)-1);

    return e0kp/(hbar*_omega);
    }
    T calce0tau(T Tl)
    {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
	
    e0kp=hbar*_omega/(exp(hbar*_omega/kb/Tl)-1)/_tau;
    //e0kp=1./(exp(hbar*_omega/kb/Tl)-1)/_tau;
    //T eref=1./(exp(hbar*_omega/kb/300)-1)/_tau;
	
    return e0kp/(hbar*_omega);
    }
    T calcde0taudT(T Tl)
    {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
	
    e0kp=(kb/_tau)*pow((hbar*_omega/kb/Tl),2)*exp(hbar*_omega/kb/Tl)/
    pow((exp(hbar*_omega/kb/Tl)-1),2);
	
    //e0kp=pow(Tl,-2)/_tau*hbar*_omega/kb*exp(hbar*_omega/kb/Tl)/
    // pow((exp(hbar*_omega/kb/Tl)-1),2);
	
    return e0kp/(hbar*_omega);
    }
    T calcde0dT(T Tl)
    {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
	
    e0kp=kb*pow((hbar*_omega/kb/Tl),2)*exp(hbar*_omega/kb/Tl)/
    pow((exp(hbar*_omega/kb/Tl)-1),2);

    //e0kp=pow(Tl,-2)*hbar*_omega/kb*exp(hbar*_omega/kb/Tl)/
    // pow((exp(hbar*_omega/kb/Tl)-1),2);
	
    return e0kp/(hbar*_omega);
    }*/
      

  //====================//
  //  Constant Cp Methods     //
  //====================//
      
  T calce0(T Tl) 
  {
    //const T hbar=6.582119e-16;  // (eV s)
    return (Tl-_Tref)*_cp;
  }
  T calce0tau(T Tl) 
  {
    //const T hbar=6.582119e-16;  // (eV s)
    return (Tl-_Tref)*_cp/_tau;
  }
  T calcde0taudT(T Tl) 
  {
    //const T hbar=6.582119e-16;  // (eV s)
    return _cp/_tau;
  }
  T calcde0dT(T Tl) 
  {
    //const T hbar=6.582119e-16;  // (eV s)
    return _cp;
  }

      
  T calcTensorPrefactor(T Tl)
  {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
	
    e0kp=exp(hbar*_omega/kb/Tl)/pow((exp(hbar*_omega/kb/Tl)-1),2)
      /_tauN/kb/Tl;
	
    return e0kp; 
  }
      
  T calcVectorPrefactor(T Tl, T e)
  {
    const T hbar=6.582119e-16;  // (eV s)
    T e0=calce0(Tl);
    e0=(e0-e)/_tauN/hbar/_omega;
    return e0;
  }
      
  T calcShifted(T Tl, T shift)
  {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K)
    T e0kp;
	
    e0kp=hbar*_omega/(exp((hbar*_omega-shift)/kb/Tl)-1);
	
    return e0kp;
  }
      
  //========================//
  //  End Constant Cp Methods     //
  //=======================//
      
  Tmode& operator=(Tmode& o)
    {
      _vg=o.getv();
      _cp=o.getcp();
      _tau=o.gettau();
      _omega=o.getomega();
      _reflections=o.getreflmap();
	  
      return *this;
    }

  T getTref() {return _Tref;}
      
  void copyPmode(Tmode& inMode)
  {
    Tvec newVg=inMode.getv();
    T newCp=inMode.getcp();
    T newTau=inMode.gettau();
    T newOmega=inMode.getomega();
    Refl_Map& newMap=inMode.getreflmap();
    int newIndex=inMode.getIndex();
    T newTref=inMode.getTref();
	
    _vg=newVg;
    _cp=newCp;
    _tau=newTau;
    _omega=newOmega;
    _reflections=newMap;
    _index=newIndex;
    _Tref=newTref;
  }
      
 private:
      
  pmode(const pmode&);
  //group velocities  (m/s)
  Tvec _vg;
      
  //specific heat
  T _cp;
      
  //relaxation rate
  T _tau;
      
  //Normal process relaxation
  T _tauN;
      
  //frequency  (rad/s)
  T _omega;
      
  //e"
  Field _efield;
      
  //shifted e" field
  Field _eShifted;
      
  //e0
  Field _e0field;
      
  //injected field
  Field _injected;
      
  //residual
  Field _residual;
      
  //residual of the injected solution, minus the injected residual
  Field _FASCorrection;
      
  //Map for specular reflections
  Refl_Map _reflections;
      
  //index for point matrix
  int _index;

  T _Tref;
};

#endif

    
