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
  
  typedef Vector<T,3> Tvec;
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
    _e0field("e0"),
    _residual("residual"),
    _reflections()
      {}
  
  Tvec getv(){return _vg;}
  T getcp() {return _cp;}
  T gettau() {return _tau;}
  T getomega() {return _omega;}
  Refl_Map& getreflmap() {return _reflections;}
  Refl_pair& getReflpair(int i) {return _reflections[i];}
  Field& getfield() {return _efield;}
  Field& gete0field() {return _e0field;}
  Field& getresid() {return _residual;}
  T calce0(T Tl)
  {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K)
    T e0kp;
    
    e0kp=hbar*_omega/(exp(hbar*_omega/kb/Tl)-1);

    return e0kp;
  }
  T calce0tau(T Tl)
  {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
    
    e0kp=hbar*_omega/(exp(hbar*_omega/kb/Tl)-1)/_tau;

    return e0kp;
  }
  T calcde0taudT(T Tl)
  {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
    
    e0kp=(kb/_tau)*pow((hbar*_omega/kb/Tl),2)*exp(hbar*_omega/kb/Tl)/
      pow((exp(hbar*_omega/kb/Tl)-1),2);

    return e0kp;
  }
  T calcde0dT(T Tl)
  {
    const T hbar=6.582119e-16;  // (eV s)
    const T kb=8.617343e-5;  // (eV/K) 
    T e0kp;
    
    e0kp=kb*pow((hbar*_omega/kb/Tl),2)*exp(hbar*_omega/kb/Tl)/
      pow((exp(hbar*_omega/kb/Tl)-1),2);
    
    return e0kp; 
  }
  
 private:
  
  pmode(const pmode&);
  //group velocities  (m/s)
  Tvec _vg;

  //specific heat
  T _cp;
  
  //relaxation rate
  T _tau;
  
  //frequency  (rad/s)
  T _omega;

  //e"
  Field _efield;

  //e0
  Field _e0field;

  //residual
  Field _residual;

  //Map for specular reflections
  Refl_Map _reflections;
};

#endif

    
