// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ARRAY_H_
#define _ARRAY_H_


#include "NumType.h"
#include "ArrayBase.h"


template <class T>
class Array : public ArrayBase
{

public:

  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  
  explicit Array(const int length) :
    ArrayBase(),
    _length(length),
    _data(new T[length]),
    _ownData(true)
  {
    logCtorVerbose("of   length %d" , _length);
  }
  
  

  Array(Array& parent, const int offset, const int length) :
    ArrayBase(),
    _length(length),
    _data(parent._data+offset),
    _ownData(false)
  {
    logCtorVerbose("with offset %d and length %d" , offset, _length);
  }

  
  virtual ~Array()
  {
    if (_ownData)
    {
        logDtorVerbose("of length %d with own data" , _length);
        delete [] _data;
    }
    else
    {
        logDtorVerbose("of length %d with offset data" , _length);
    }
  }

  void resize(const int newLength)
  {
    if (_ownData)
    {
        T* newData = new T[newLength];
        if (_data)
        {
	  int len;
	  if(newLength<_length)
	    len=newLength;
	  else
	    len=_length;

	  for(int i=0; i<len; i++)
	    newData[i] = _data[i];
	  delete [] _data;
        }
        _data = newData;
        _length = newLength;
    }
    else
      throw CException("cannot resize offset array");
  }
  
  DEFINE_TYPENAME("Array<" + NumTypeTraits<T>::getTypeName() + ">");
  
  virtual int getDimension() const
  {
    return NumTypeTraits<T>::getDimension() + 1;
  };

  int getLength() const {return _length;}

  virtual PrimType getPrimType() const
  {
    return NumTypeTraits<T>::getPrimType();
  }
  
  virtual void getShape(int* shp) const
  {
    *shp = _length;
    NumTypeTraits<T>::getShape(shp+1);
  }

  
  T& operator[](int n) {return _data[n];}
  const T& operator[](int n) const {return _data[n];}

  void operator=(const T& x)
  {
    for(int i=0;i<_length;i++)
      _data[i] = x;
  }
  
  void operator=(const Array& o)
  {
    for(int i=0;i<_length;i++)
      _data[i] = o._data[i];
  }

  
  virtual ArrayBase& operator+=(const ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (o._length == _length)
    {
        for(int i=0;i<_length;i++)
          _data[i] += o._data[i];
    }
    else
      throw CException("invalid array for operator +=");          
    return *this;
  }

  virtual ArrayBase& operator-=(const ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (o._length == _length)
    {
        for(int i=0;i<_length;i++)
          _data[i] -= o._data[i];
    }
    else
      throw CException("invalid array for operator -=");          
    return *this;
  }
  
  virtual ArrayBase& operator/=(const ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (o._length == 1)
    {
        for(int i=0;i<_length;i++)
          _data[i] /= o._data[0];
    }
    else if (o._length == _length)
    {
        for(int i=0;i<_length;i++)
          _data[i] /= o._data[i];
    }
    else
      throw CException("invalid array for operator /=");          
    return *this;
  }

  virtual void setMax(const  ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (_length == 1 && o._length==1)
    {
        for(int i=0;i<_length;i++)
          NumTypeTraits<T>::setMax(_data[i], o._data[i]);
    }
    else
      throw CException("invalid array for setMax");          
  }

  virtual ArrayBase& safeDivide(const  ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (_length == 1 && o._length==1)
    {
        for(int i=0;i<_length;i++)
          NumTypeTraits<T>::safeDivide(_data[i], o._data[i]);
    }
    else
      throw CException("invalid array for safeDivide");          
    return *this;
  }

  virtual ArrayBase& normalize(const  ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (_length == 1 && o._length==1)
    {
        for(int i=0;i<_length;i++)
          NumTypeTraits<T>::normalize(_data[i], o._data[i]);
    }
    else
      throw CException("invalid array for normalize");          
    return *this;
  }


  virtual ArrayBase& operator*=(const ArrayBase& obase)
  {
    const Array& o = dynamic_cast<const Array& >(obase);
    if (o._length == 1)
    {
        for(int i=0;i<_length;i++)
          _data[i] *= o._data[0];
    }
    else if (o._length == _length)
    {
        for(int i=0;i<_length;i++)
          _data[i] *= o._data[i];
    }
    else
      throw CException("invalid array for operator *=");          
    return *this;
  }

  virtual bool operator<(const double tolerance) const
  {
    if (_length == 1)
    {
        return (_data[0] <tolerance);
    }
    else
      throw CException("invalid array for operator<");          
  }

  virtual void print(ostream& os) const 
  {
    if (_length > 0)
    {
       for ( int i = 0; i < _length; i++ )
         if ( i < _length-1)
            os <<_data[i] << endl;
         else 
            os <<_data[i];
    }
    else
      throw CException("invalid array for print");          
  }

  // this += alpha*x
  virtual ArrayBase& saxpy(const ArrayBase& alphabase, const ArrayBase& xbase)
  {
    const Array& alpha = dynamic_cast<const Array& >(alphabase);
    const Array& x = dynamic_cast<const Array& >(xbase);
    
    if (alpha._length == 1 && x._length == _length)
    {
        for(int i=0;i<_length;i++)
          _data[i] += alpha._data[0]*x._data[i];
    }
    else
      throw CException("invalid arrays for saxpy");          
    return *this;
  }

  
  // this -= alpha*x; defined to avoid having to create a negative of array
  virtual ArrayBase& msaxpy(const ArrayBase& alphabase, const ArrayBase& xbase)
  {
    const Array& alpha = dynamic_cast<const Array& >(alphabase);
    const Array& x = dynamic_cast<const Array& >(xbase);
    
    if (alpha._length == 1 && x._length == _length)
    {
        for(int i=0;i<_length;i++)
          _data[i] -= alpha._data[0]*x._data[i];
    }
    else
      throw CException("invalid arrays for msaxpy");          
    return *this;
  }

  virtual void* getData() const {return  _data;}
  virtual int getDataSize() const
  {
    return  _length*NumTypeTraits<T>::getDataSize();
  }

  virtual void zero()
  {
    memset(_data,0,getDataSize());
  }
  
  void disownData() const {_ownData = false;}

  

  virtual shared_ptr<ArrayBase>
  dotWith(const ArrayBase& abase, const int lengthToUse) const
  {
    const Array& a = dynamic_cast<const Array&>(abase);
    T sum(NumTypeTraits<T>::getZero());
    for(int i=0; i<lengthToUse; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum,_data[i],a._data[i]);
    shared_ptr<Array> nPtr(new Array(1));
    (*nPtr)[0] = sum;
    return nPtr;
  }

  virtual shared_ptr<ArrayBase>
  getOneNorm(const int lengthToUse) const
  {
    T sum(NumTypeTraits<T>::getZero());
    for(int i=0; i<lengthToUse; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum,_data[i]);
    shared_ptr<Array> nPtr(new Array(1));
    (*nPtr)[0] = sum;
    return nPtr;
  }

  virtual shared_ptr<ArrayBase>
  reduceSum() const
  {
    T_Scalar sum(NumTypeTraits<T_Scalar>::getZero());
    for(int i=0; i<_length; i++)
      NumTypeTraits<T>::reduceSum(sum,_data[i]);
    shared_ptr<Array<T_Scalar> > nPtr(new Array<T_Scalar>(1));
    (*nPtr)[0] = sum;
    return nPtr;
  }

  virtual void
  setSum(const ArrayBase& sumBase)
  {
    const Array<T_Scalar>& sum =
      dynamic_cast<const Array<T_Scalar>& >(sumBase);
    if (sum.getLength() == 1)
    {
        const T_Scalar& s = sum[0];
        for(int i=0; i<_length; i++)
          _data[i]=s;
    }
  }

  virtual void
  scatter(ArrayBase& other_, const ArrayBase& indices_, const int offset=0) const
  {
    const Array<int>& indices = dynamic_cast<const Array<int>& >(indices_);
    Array& other = dynamic_cast<Array& >(other_);
    for(int ii=0; ii<indices.getLength(); ii++)
        other[offset+ii] = _data[indices[ii]];
  }

  virtual void
  gather(const ArrayBase& other_, const ArrayBase& indices_, const int offset=0)
  {
    const Array<int>& indices = dynamic_cast<const Array<int>& >(indices_);
    const Array& other = dynamic_cast<const Array& >(other_);
    for(int ii=0; ii<indices.getLength(); ii++)
      _data[indices[ii]] = other[offset+ii];
  }

  virtual shared_ptr<Array>
  getSubset(const Array<int>& indices)
  {
    shared_ptr<Array> subPtr(new Array(indices.getLength()));
    Array& sub = *subPtr;

    for(int ii=0; ii<indices.getLength(); ii++)
    {
        sub[ii] = _data[indices[ii]];
    }

    return subPtr;
  }

  virtual void
  setSubsetFromSubset(const ArrayBase& other_, const ArrayBase& fromIndices_,
                      const ArrayBase& toIndices_)
  {
    const Array<int>& toIndices = dynamic_cast<const Array<int>& >(toIndices_);
    const Array<int>& fromIndices = dynamic_cast<const Array<int>& >(fromIndices_);
    const Array& other = dynamic_cast<const Array& >(other_);
    for(int ii=0; ii<fromIndices.getLength(); ii++)
    {
      _data[toIndices[ii]] = other[fromIndices[ii]];
    }
  }

  virtual void
  copyFrom(const IContainer& oc)
  {
    const Array& other = dynamic_cast<const Array&>(oc);
    operator=(other);
  }

  virtual void
  copyPartial(const IContainer& oc, const int iBeg, const int iEnd)
  {
    const Array& other = dynamic_cast<const Array&>(oc);
    for(int i=iBeg;i<iEnd;i++)
      _data[i] = other._data[i];
  }

  virtual void
  zeroPartial(const int iBeg, const int iEnd)
  {
    for(int i=iBeg;i<iEnd;i++)
      _data[i] = NumTypeTraits<T>::getZero();
  }

  virtual void limit(const double min, const double max)
  {
    if (_length == 1)
    {
        ArrayScalarTraits<T>::limit(_data[0],min, max);
    }
    else
      throw CException("invalid array for limit");
  }

  virtual shared_ptr<ArrayBase>
  operator-() const
  {
    if (_length == 1)
    {
        Array* nPtr(new Array(1));
        (*nPtr)[0] = -_data[0];
        return shared_ptr<ArrayBase>(nPtr);
    }
    throw CException("invalid  array for operator-");
  }


  virtual void
  inject(IContainer& coarseI, const IContainer& coarseIndexI, const int length) const
  {
    Array& coarse = dynamic_cast<Array&>(coarseI);
    const Array<int>& coarseIndex = dynamic_cast<const Array<int>& >(coarseIndexI);

    for(int i=0; i<length; i++)
      if (coarseIndex[i]>=0)
        coarse[coarseIndex[i]] += _data[i];
  }

  virtual void
  correct(const IContainer& coarseI, const IContainer& coarseIndexI,
          const IContainer* scaleI, const int length)
  {
    const Array& coarse = dynamic_cast<const Array&>(coarseI);
    const Array<int>& coarseIndex = dynamic_cast<const Array<int>& >(coarseIndexI);

    if (scaleI)
    {
        const Array& scaleArray =
          dynamic_cast<const Array&>(*scaleI);
        if (scaleArray.getLength() == 1)
        {
            const T& scale = scaleArray[0];
            //cout << "correction scale" << scale << endl;

            for(int i=0; i<length; i++)
              _data[i] += coarse[coarseIndex[i]]*scale;

        }
        else
          throw CException("invalid scale array");
    }
    else
    {
        for(int i=0; i<length; i++)
          if (coarseIndex[i]>=0)
            _data[i] += coarse[coarseIndex[i]];
    }
  }

  virtual shared_ptr<IContainer>
  newClone() const
  {
    return shared_ptr<Array>(new Array(_length));
  }

  virtual shared_ptr<ArrayBase>
  newSizedClone(const int size) const
  {
    return shared_ptr<Array>(new Array(size));
  }
    

  virtual shared_ptr<IContainer>
  newCopy() const
  {
    shared_ptr<Array> c(new Array(_length));
    *c = *this;
    return c;
  }
  
  virtual shared_ptr<ArrayBase>
  createOffsetArray(const int offset, const int length)
  {
    return shared_ptr<Array>(new Array(*this,offset,length));
  }

private:
  Array(const Array&);

  int _length;
  T* _data;
  mutable bool _ownData;
};



template<class T>
class ConstAsArray
{
public:
  ConstAsArray(const T& c) :
    _c(c)
  {}

  const T& operator[](int n) const {return _c;}

private:
  const T& _c;
};


template<class T>
shared_ptr<Array<T> >
arrayFromVector(const vector<T>& v)
{
  const int count = v.size();
  Array<T>* a = new Array<T>(count);
  for(int i=0; i<count; i++)
    (*a)[i] = v[i];
  return shared_ptr<Array<T> >(a);
}

#endif
