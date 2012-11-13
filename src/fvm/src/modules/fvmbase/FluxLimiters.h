// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLUXLIMITERS_H_
#define _FLUXLIMITERS_H_

#include <math.h>

struct minModLim
{

  template<class T>
  T operator()(const T r)
  {return max(0., min(1., r));}

};

struct vanLeer
{
  template<class T>
  T operator()(const T a, const T b) const
  {
    if (b!=0.)
      return (a/b+fabs(a/b))/(1+fabs(a/b));
    else 
      return 0.;
  }
  
  template<class T>
  T operator()(const T r) const
  {
    if (r!=0.)
      return (r+fabs(r))/(1+fabs(r));
    else 
      return 0.;
  }

};

struct superbee
{
  template<class T>
  T operator()(const T a, const T b) const
  {
    if(b!=0.)
      return max(0., max(min(1., 2.*a/b), min(2., a/b)));
  }
};

struct vanAlbada1
{
  template<class T>
  T operator()(const T r)
  {return (pow(r,2)+r)/(pow(r,2)+1.);}
};

struct ospre
{

  template<class T>
  T operator()(const T a, const T b) const
  {
    if(a>0 && b!=0.)
      return 1.5*(pow(a/b,2)+a/b)/(pow(a/b,2)+a/b+1.);
    return 0.;
  }

};

template<class T, class LimitFunc>
  void computeLimitCoeff(T& lc, const T x, const T& dx, 
			 const T& min, const T& max, const LimitFunc& f)
{
  T thislc = NumTypeTraits<T>::getUnity();
  if (dx>0)
    thislc = f(dx,max-x);
  else
    thislc = f(-dx,x-min);
  if (thislc < lc)
    lc = thislc;
};

template<class T, class LimitFunc>
  void computeLimitCoeff2(T& lc, const T x, const T& dx, 
			 const T& min, const T& max, const LimitFunc& f)
{
  T thislc = NumTypeTraits<T>::getUnity();
  if (dx>0)
    thislc = f(dx/(max-x));
  else
    thislc = f(-dx/(x-min));
  if (thislc < lc)
    lc = thislc;
};


#endif
