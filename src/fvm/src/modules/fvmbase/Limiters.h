// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LIMITERS_H_
#define _LIMITERS_H_

#include <math.h>

struct DefaultLimiter
{
  template<class T>
  T operator()(const T& a, const T& b) const
  {
    /*
    if (b < T(2.0)*a)
      return (NumTypeTraits<T>::getUnity());

    const T eps(1e-16);
    const T c = b*b + T(4.)*a*a;
    return ((b*c+eps) / (b*b*b + a*c + eps));
    */
  }
};

struct SuperbeeLimiter
{
 template<class T>
 T operator()(const T a, const T b) const
  {
    if(b!=0.)
      {
	//T r = (4.*(a/b))-1.;
	//return max(0., max(min(1., 2.*a/b), min(2., a/b)));
	//return max(0., max(min(1., a/b), min(2., a/(2.*b))));
	//return max(0., max(min(1., 2.*r), min(2., r)));
	//if (b < T(2.0)*a)
	//return (NumTypeTraits<T>::getUnity());

	const T eps(1e-16);
	const T c = b*b + T(4.)*a*a;
	return ((b*c+eps) / (b*b*b + a*c + eps));
      }
  }
};

/*
template<class T, class LimitFunc>
  void computeLimitCoeff(T& lc, const T x,
			 const T& dx, const T& min, const T& max, const LimitFunc& f)
{
  T thislc = NumTypeTraits<T>::getUnity();
  if (dx>0)
    thislc = f(dx,max-x);
  else
    thislc = f(-dx,x-min);
  if (thislc < lc)
    lc = thislc;
}
*/



#endif
