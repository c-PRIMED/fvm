#ifndef _KVOL_H_
#define _KVOL_H_

#include <vector>
#include "Vector.h"
#include "Field.h"
#include "pmode.h"

template<class T>
class kvol
{
  
 public:

  typedef Vector<T,3> Tvec;
  typedef pmode<T> Tmode;
  typedef shared_ptr<Tmode> Tmodeptr;
  typedef vector<Tmodeptr> Modes;

 kvol(Tmodeptr mode,T dk3):     //used in gray approx.
  _dk3(dk3),
    _modenum(1),
    _modes(1,mode)
      {}
  
  Tvec getkvec() {return _Kvector;}
  void setkvec(Tvec K) {_Kvector=K;}
  void setdk3(T dk3) {_dk3=dk3;}
  T getdk3() {return _dk3;}
  int getmodenum() {return _modenum;}
  Tmode& getmode(int n) const {return *_modes[n];}

 private:

  kvol(const kvol&);

  //weight factor (1/m3)
  T _dk3;

  //K vector
  Tvec _Kvector;

  //number of modes
  int _modenum;

  //phonon modes
  Modes _modes;
  
};

#endif
