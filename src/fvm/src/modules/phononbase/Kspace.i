%{
#include "Kspace.h"
  %}


template<class T>
class Kspace
{

 public:

  typedef Vector<T,3> Tvec;
  typedef pmode<T> Tmode;
  typedef shared_ptr<Tmode> Tmodeptr;
  typedef kvol<T> Tkvol;
  typedef shared_ptr<Tkvol> Kvolptr;
  typedef vector<Kvolptr> Volvec;
  typedef typename Tmode::Reflection Reflection;
  typedef typename Tmode::Reflptr Reflptr;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef typename Tmode::Refl_Map Refl_Map;

  Kspace(T a, T tau, T vgmag, T omega, int ntheta, int nphi) ;
  int getlength();
  int gettotmodes();
  T getDK3() ;
  T calcSpecificHeat(T Tl);
  
 private:
  
  Kspace(const Kspace&);
  //num volumes
  int _length;
  Volvec _Kmesh;
  T _totvol;    //total Kspace volume
  
};


%template(KspaceA) Kspace< ATYPE_STR >;
