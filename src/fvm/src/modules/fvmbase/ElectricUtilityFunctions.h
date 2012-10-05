// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ELECTRICUTILITYFUNCTIONS_H_
#define _ELECTRICUTILITYFUNCTIONS_H_

/*** a collection of utility functions that are used in dielectric charging model ***/

#include "PhysicsConstant.h"
#include "Array.h"
#include "Vector.h"
#include "CException.h"


typedef Vector<double, 3> VectorD3;

template<class T>
const T FermiFunction ( const T& energy, const T& fermilevel, const T& temperature)
{
  return 1. / (1. + exp(QE * (energy - fermilevel) / (K_SI * temperature) ));
}

template<class T>
const T  ElectronSupplyFunction (const T& energy, const T& fermilevel, const T& temperature)
{
  T exponent, power, supply;
    
  power = -QE * (energy - fermilevel) / (K_SI * temperature);
  exponent = exp(power);
    
  if (exponent <= 0.01)
    supply = K_SI * temperature * (exponent - pow(exponent, 2.0)/2.0 + pow(exponent, 3.0)/3.0 - pow(exponent, 4.0)/4.0);
    
  if (power >= 10.0)
    supply = K_SI * temperature * power;
  
  else
    supply = K_SI * temperature *log(1+exponent);
	

  //supply = K_SI * temperature *log(1+exponent);
  return supply;
}




template<class T>
const T PositiveValueOf(T input) {
  if (input>=0)
    return input;
  else
    return T(0);
}

template<class T>
const T getMembraneVoltage ( const T& currentTime)
{
  T volt;
  volt = 1.0;
  return volt;
}

double SignOf(double x) {
  if (x>0) return 1;
  else if (x<0) return -1;
  else return 0; 
}

double findMin(double x1, double x2, double x3, double x4)
{
  double min = 1000000000;
  if (x1 <= min)
    min = x1;
  if (x2 <= min)
    min = x2;
  if (x3 <= min)
    min = x3;
  if (x4 <= min)
    min = x4;
  return min;
}
  
double findMax(double x1, double x2, double x3, double x4)
{
  double max = -1000000000;
  if (x1 >= max)
    max = x1;
  if (x2 >= max)
    max = x2;
  if (x3 >= max)
    max = x3;
  if (x4 >= max)
    max = x4;
  return max;
}

double distanceFromPointToLine(VectorD3 p1, VectorD3 p2, VectorD3 M)
{
  // 3D point p1 and p2 forms a line L
  // find the shortest distance from point M to line L
  // assume point H is the projection of M on line L
  // so the shortest distance is the length of MH
  
  VectorD3 v = p2 - p1;
  if (mag(v) == 0) 
    throw CException ("column center line start and end points are the same!!!");
  
  v /= mag(p2-p1);
  VectorD3 PM = p1 - M;
  VectorD3 PH = dot(PM, v) * v;
  VectorD3 HM = PM - PH;
  double distance = mag(HM);
  return distance; 
}

VectorD3 projectionFromPointToLine(VectorD3 p1, VectorD3 p2, VectorD3 M)
{
  VectorD3 v = p2 - p1;
  if (mag(v) == 0) 
    throw CException ("column center line start and end points are the same!!!");
  
  v /= mag(p2-p1);
  VectorD3 PM = p1 - M;
  VectorD3 PH = dot(PM, v) * v;
  VectorD3 HM = PM - PH;
  VectorD3 H = HM + M;
  return H;
}
#endif
