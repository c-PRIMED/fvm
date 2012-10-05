// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _TRACTIONVAL_H_
#define _TRACTIONVAL_H_

#include <map>
#include "misc.h"
#include "GeomFields.h"
#include "Vector.h" 
#include "NumType.h"
#include "StorageSite.h"
#include "FloatVarDict.h"
class Field;

/*
 * Nine stress components should be provided as arguments in the order:
 * TractionXX,TractionXY,TractionXZ,TractionYX,TractionYY,TractionYZ,
 * TractionZX,TractionZY,TractionZZ. These components can be provided either
 * as fields or as values. However, all of TractionXX,TractionXY,and TractionXZ
 * should be either fields or values. Similarly, all of TractionYX,TractionYY,
 * and TractionYZ should either be fields or values. The same holds true for
 * TractionZX,TractionZY,and TractionZZ.
 */

template <class T>
class TractionValEvaluator
{
public:
  TractionValEvaluator(const FloatVal<T>& fval0,
                    const FloatVal<T>& fval1,
                    const FloatVal<T>& fval2,
		    const FloatVal<T>& fval3,
		    const FloatVal<T>& fval4,
		    const FloatVal<T>& fval5,
		    const FloatVal<T>& fval6,
		    const FloatVal<T>& fval7,
		    const FloatVal<T>& fval8,
		    const GeomFields& geomFields,
                    const StorageSite& site) :
    _faceArea(dynamic_cast<const Array<Vector<T,3> >&>(geomFields.area[site])),
    _faceAreaMag(dynamic_cast<const Array<T>&>(geomFields.areaMag[site])),
    _isField0(fval0.field != 0),
    _isField1(fval4.field != 0),
    _isField2(fval8.field != 0),
    _arrayPtr0(),
    _arrayPtr1(),
    _arrayPtr2(),
    _arrayPtr3(),
    _arrayPtr4(),
    _arrayPtr5(),
    _arrayPtr6(),
    _arrayPtr7(),
    _arrayPtr8()
  {
    if (_isField0)
    {
        _arrayPtr0 = dynamic_pointer_cast<Array<T> >(fval0.field->getArrayPtr(site));
	_arrayPtr1 = dynamic_pointer_cast<Array<T> >(fval1.field->getArrayPtr(site));
	_arrayPtr2 = dynamic_pointer_cast<Array<T> >(fval2.field->getArrayPtr(site));
    }
    else
    {
        _constant0[0] = fval0.constant;
        _constant0[1] = fval1.constant;
        _constant0[2] = fval2.constant;
    }
    if (_isField1)
    {
        _arrayPtr3 = dynamic_pointer_cast<Array<T> >(fval3.field->getArrayPtr(site));
        _arrayPtr4 = dynamic_pointer_cast<Array<T> >(fval4.field->getArrayPtr(site));
        _arrayPtr5 = dynamic_pointer_cast<Array<T> >(fval5.field->getArrayPtr(site));
    }
    else
    {
        _constant1[0] = fval3.constant;
        _constant1[1] = fval4.constant;
        _constant1[2] = fval5.constant;
    }      
    if (_isField2)
    {
        _arrayPtr6 = dynamic_pointer_cast<Array<T> >(fval6.field->getArrayPtr(site));
	_arrayPtr7 = dynamic_pointer_cast<Array<T> >(fval7.field->getArrayPtr(site));
	_arrayPtr8 = dynamic_pointer_cast<Array<T> >(fval8.field->getArrayPtr(site));
    }
    else
    {
        _constant2[0] = fval6.constant;
        _constant2[1] = fval7.constant;
        _constant2[2] = fval8.constant;
    }
  }

  const Vector<T,3>& operator[](const int i) const
  {
   
    Vector<T,3> en(_faceArea[i]/_faceAreaMag[i]);
    Vector<T,3> traction0(NumTypeTraits<Vector<T,3> >::getZero());
    Vector<T,3> traction1(NumTypeTraits<Vector<T,3> >::getZero());
    Vector<T,3> traction2(NumTypeTraits<Vector<T,3> >::getZero());

    if (_isField0)
    {
        traction0[0] = (*_arrayPtr0)[i];
	traction0[1] = (*_arrayPtr1)[i];
	traction0[2] = (*_arrayPtr2)[i];
	_constant[0] = dot(traction0,en);
    }
    else
    {
        _constant[0] = dot(_constant0,en);
    }
    if (_isField1)
    {
        traction1[0] = (*_arrayPtr3)[i];
	traction1[1] = (*_arrayPtr4)[i];
	traction1[2] = (*_arrayPtr5)[i];
	_constant[1] = dot(traction1,en);
    }
    else
    {
        _constant[1] = dot(_constant1,en);
    }
    if (_isField2)
    {
        traction2[0] = (*_arrayPtr6)[i];
	traction2[1] = (*_arrayPtr7)[i];
	traction2[2] = (*_arrayPtr8)[i];
	_constant[2] = dot(traction2,en);
    }
    else
    {
        _constant[2] = dot(_constant2,en);
    }

    return _constant;
  }
  
private:
  mutable Vector<T,3> _constant0;
  mutable Vector<T,3> _constant1;
  mutable Vector<T,3> _constant2;
  mutable Vector<T,3> _constant;
  const Array<Vector<T,3> >& _faceArea;
  const Array<T>& _faceAreaMag; 
  const bool _isField0;
  const bool _isField1;
  const bool _isField2;
  shared_ptr<Array<T> > _arrayPtr0;
  shared_ptr<Array<T> > _arrayPtr1;
  shared_ptr<Array<T> > _arrayPtr2;
  shared_ptr<Array<T> > _arrayPtr3;
  shared_ptr<Array<T> > _arrayPtr4;
  shared_ptr<Array<T> > _arrayPtr5;
  shared_ptr<Array<T> > _arrayPtr6;
  shared_ptr<Array<T> > _arrayPtr7;
  shared_ptr<Array<T> > _arrayPtr8;
};

#endif
