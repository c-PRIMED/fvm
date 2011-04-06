#ifndef _KSPACE_H_
#define _KSPACE_H_

#include <vector>
#include <math.h>
#include "Vector.h"
#include "pmode.h"
#include "kvol.h"

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

 Kspace(T tau,T vgmag,T cp,int ntheta,int nphi) :
  _length(ntheta*nphi),
    _Kmesh(),
    _totvol(0.),
    _Covertau(0.)
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
		Tmodeptr modeptr=shared_ptr<Tmode>(new Tmode(vg,cp,tau));
		Kvolptr volptr=shared_ptr<Tkvol>(new Tkvol(modeptr,dk3));
		_Kmesh.push_back(volptr);
		_Covertau+=cp/tau*dk3;
		_totvol+=dk3;
	      }
	  }
	_Covertau=_Covertau/_totvol;
      }
  
  //void setvol(int n,Tkvol k) {*_Kmesh[n]=k;}
  Tkvol& getkvol(int n) const {return *_Kmesh[n];}
  int getlength() const {return _length;}
  int gettotmodes()
  {
    return (_Kmesh[0]->getmodenum())*_length;
  }
  T getDK3() const {return _totvol;}
  T getCovertau() const {return _Covertau;}

 private:

  Kspace(const Kspace&);
  //num volumes
  int _length;
  Volvec _Kmesh;
  T _totvol;    //total Kspace volume
  T _Covertau;
  

};


#endif
