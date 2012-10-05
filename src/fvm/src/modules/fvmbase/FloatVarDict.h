// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _FLOATVARDICT_H_
#define _FLOATVARDICT_H_

#include <map>
#include "misc.h"
class Field;

/**
 * FloatVal allows specification of any floating point input as either
 * a single number or as a Field. When a Field is specified it must
 * contain Array<T> for any storage site for which the input needs a
 * value.
 * 
 */

template<class T>
struct FloatVal
{
  FloatVal(T constant_) :
    constant(constant_),
    field(0)
  {}

  FloatVal(Field* field_) :
    constant(0),
    field(field_)
  {}

  T constant;
  Field *field;
};

/**
 * A dictionary for specifying floating point inputs -- used for
 * boundary conditions as well as model inputs.  Derived classes
 * define the keys and default values.
 * 
 */

template<typename T>
class FloatVarDict : public map<string, FloatVal<T> >
{
public:

  typedef map<string,FloatVal<T> > T_Parent;

  bool hasVar(const string varName) const
  {
    typename T_Parent::const_iterator pos = this->find(varName);
    return pos != this->end();
  }
  
  T operator[](const string varName) const
  {
    typename T_Parent::const_iterator pos = this->find(varName);
    if (pos != this->end())
    {
        if (pos->second.field)
          throw CException(varName + " value is a Field");
        return pos->second.constant;
    }
    throw CException("uknown var " + varName);
  }

  bool isField(const string varName) const
  {
    typename T_Parent::const_iterator pos = this->find(varName);
    if (pos != this->end())
      return (pos->second.field!=0);
    throw CException("uknown var " + varName);
  }
  
  Field& getField(const string varName) const
  {
    typename T_Parent::const_iterator pos = this->find(varName);
    if (pos != this->end())
      return *(pos->second.field);
    throw CException("uknown var " + varName);
  }
  
  FloatVal<T> getVal(const string varName) const
  {
    typename T_Parent::const_iterator pos = this->find(varName);
    if (pos != this->end())
      return pos->second;
    throw CException("uknown var " + varName);
  }
  
  protected:
  void defineVar(const string varName, const T defaultValue)
  {
    this->insert(make_pair(varName,defaultValue));
  }
};

/**
 * In dealing with inputs defined by FloatingVarDict we can use
 * FloatValEvaluator as an array to read values at each index in the
 * specified site. When the input value is a field this will return
 * the corresponding array value from the Field; otherwise the
 * constant value is returned for all indices.
 * 
 */

template <class T>
class FloatValEvaluator
{
public:
  FloatValEvaluator(const FloatVal<T>& fval, const StorageSite& site) :
    _constant(fval.constant),
    _isField(fval.field != 0),
    _arrayPtr()
  {
    if (_isField)
      _arrayPtr = dynamic_pointer_cast<Array<T> >(fval.field->getArrayPtr(site));
    
  }

  FloatValEvaluator(const T val) :
    _constant(val),
    _isField(false),
    _arrayPtr(0)
  {}

  T operator[](const int i) const
  {
    return (_isField ? (*_arrayPtr)[i] : _constant);
  }
  
private:
  const T _constant;
  const bool _isField;
  shared_ptr<Array<T> > _arrayPtr;
};


/**
 * Since all the values handled by FloatingVarDict are single floating
 * point numbers, we need this template specialization of
 * FloatValEvaluator for cases where we are expecting a Vector<T,3>
 * variable. It always returns a vector whose components may be either
 * the constant value or value from a Field array. This allows us to
 * specify a Field only for say the y component of a vector input and
 * have constants for x and z components.
 * 
 */

template <class T>
class FloatValEvaluator<Vector<T,3> >
{
public:
  FloatValEvaluator(const FloatVal<T>& fval0,
                    const FloatVal<T>& fval1,
                    const FloatVal<T>& fval2,
                    const StorageSite& site) :
    _isField0(fval0.field != 0),
    _isField1(fval1.field != 0),
    _isField2(fval2.field != 0),
    _arrayPtr0(),
    _arrayPtr1(),
    _arrayPtr2()
  {
    if (_isField0)
      _arrayPtr0 = dynamic_pointer_cast<Array<T> >(fval0.field->getArrayPtr(site));
    else
      _constant[0] = fval0.constant;

    if (_isField1)
      _arrayPtr1 = dynamic_pointer_cast<Array<T> >(fval1.field->getArrayPtr(site));
    else
      _constant[1] = fval1.constant;
      
    if (_isField2)
      _arrayPtr2 = dynamic_pointer_cast<Array<T> >(fval2.field->getArrayPtr(site));
    else
      _constant[2] = fval2.constant;
  }

  FloatValEvaluator(const Vector<T,3>& val) :
    _constant(val),
    _isField0(false),
    _isField1(false),
    _isField2(false),
    _arrayPtr0(0),
    _arrayPtr1(0),
    _arrayPtr2(0)
  {}

  const Vector<T,3>& operator[](const int i) const
  {
    if (_isField0) _constant[0] = (*_arrayPtr0)[i];
    if (_isField1) _constant[1] = (*_arrayPtr1)[i];
    if (_isField2) _constant[2] = (*_arrayPtr2)[i];
    
    return _constant;
  }
  
private:
  mutable Vector<T,3> _constant;
  const bool _isField0;
  const bool _isField1;
  const bool _isField2;
  shared_ptr<Array<T> > _arrayPtr0;
  shared_ptr<Array<T> > _arrayPtr1;
  shared_ptr<Array<T> > _arrayPtr2;
};

#endif
