// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLUXLIMITERS_H_
#define _FLUXLIMITERS_H_

#include <math.h>

template<class T>
T minMod(const T r)
{return max(0., min(1., r));}

template<class T>
T vanLeer(const T r)
{return (r+fabs(r))/(1+fabs(r));}

template<class T>
T superbee(const T r)
{return max(0., max(min(1., 2.*r), min(2., r)));}

template<class T>
T vanAlbada1(const T r)
{return (pow(r,2)+r)/(pow(r,2)+1.);}

template<class T>
T ospre(const T r)
{
  if(r>0)
    return 1.5*(pow(r,2)+r)/(pow(r,2)+r+1.);
  return 0.;
}


#endif
