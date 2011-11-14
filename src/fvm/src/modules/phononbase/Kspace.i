%{
#include "Kspace.h"
  %}

%include "std_vector.i"

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
  typedef Kspace<T> TKspace;
  typedef typename Tmode::Reflection Reflection;
  typedef typename Tmode::Reflptr Reflptr;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef typename Tmode::Refl_Map Refl_Map;

  Kspace(T a, T tau, T vgmag, T omega, int ntheta, int nphi);
  Kspace(const char* filename,const int dimension);
  int getlength();
  int gettotmodes();
  T getDK3() ;
  T calcSpecificHeat(T Tl);

  %extend{
    Kspace<T>* getSelfPointer() {return self;}
    std::vector<Kspace<T>*> MakeList()
      {
	std::vector<Kspace<T>*> newList;
	newList.push_back(self);
	return newList;
      }
  }
  
 private:
  
  Kspace(const Kspace&);
  //num volumes
  int _length;
  Volvec _Kmesh;
  T _totvol;    //total Kspace volume
  
};


%template(KspaceA) Kspace< ATYPE_STR >;

typedef std::vector<Kspace<ATYPE_STR>*> KList;
%template(KList) std::vector<Kspace<ATYPE_STR>*>;
