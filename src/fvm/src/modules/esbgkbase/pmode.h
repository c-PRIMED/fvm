#ifndef _PMODE_H_
#define _PMODE_H_

#include "Vector.h"
#include "Field.h"
#include <string>

using namespace std;

template<class T>
class pmode
{
  
 public:
  
  typedef Vector<T,3> Tvec;
  
 pmode(Tvec vg,T cp,T tau):
  _vg(vg),
    _cp(cp),
    _tau(tau),
    _efield("edoubleprime")
      {}
  
  Tvec getv(){return _vg;}
  T getcp() {return _cp;}
  T gettau() {return _tau;}
  T getomega() {return _omega;}
  Field& getfield() {return _efield;}
  
 private:
  
  //group velocities  (m/s)
  Tvec _vg;

  //specific heat
  T _cp;
  
  //relaxation rate
  T _tau;
  
  //frequency  (rad/s)
  T _omega;

  //field for e"
  Field _efield;
  
};

#endif

    
