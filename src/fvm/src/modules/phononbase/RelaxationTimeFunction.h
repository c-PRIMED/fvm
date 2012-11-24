#ifndef _RELAXATIONTIMEFUNCTION_H_
#define _RELAXATIONTIMEFUNCTION_H_

#include <math.h>

template<class T>
struct RelTimeFun
{

  RelTimeFun()
  {
    this->_A=1.;
    this->_B=1.;
    this->_C=1.;
    this->_constTau=true;
  }

  RelTimeFun(const T A, const T B, const T C)
  {
    this->_A=A;
    this->_B=B;
    this->_C=C;
  }

  void update(const T w, const T Tl, T& tau)
  {
    if(!_constTau)
      {
	T umk=_B*Tl*pow(w,2.)*exp(-_C/Tl);
	T imp=_A*pow(w,4.);
	tau=1./(umk+imp);
      }
  }

  T _A;
  T _B;
  T _C;
  bool _constTau;
};



#endif
