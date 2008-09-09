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
  
  Array(const int length) :
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
            for(int i=0; i<_length; i++)
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
  dotWith(const ArrayBase& abase) const
  {
    const Array& a = dynamic_cast<const Array&>(abase);
    T sum(NumTypeTraits<T>::getZero());
    for(int i=0; i<_length; i++)
      NumTypeTraits<T>::accumulateDotProduct(sum,_data[i],a._data[i]);
    shared_ptr<Array> nPtr(new Array(1));
    (*nPtr)[0] = sum;
    return nPtr;
  }

  virtual shared_ptr<ArrayBase>
  getOneNorm() const
  {
    T sum(NumTypeTraits<T>::getZero());
    for(int i=0; i<_length; i++)
      NumTypeTraits<T>::accumulateOneNorm(sum,_data[i]);
    shared_ptr<Array> nPtr(new Array(1));
    (*nPtr)[0] = sum;
    return nPtr;
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

  virtual IContainer&
  operator=(const IContainer& oc)
  {
    const Array& other = dynamic_cast<const Array&>(oc);
    operator=(other);
    return *this;
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

  virtual void
  inject(IContainer& coarseI, const IContainer& coarseIndexI, const int length) const
  {
    Array& coarse = dynamic_cast<Array&>(coarseI);
    const Array<int>& coarseIndex = dynamic_cast<const Array<int>& >(coarseIndexI);

    for(int i=0; i<length; i++)
      coarse[coarseIndex[i]] += _data[i];
  }

  virtual void
  correct(const IContainer& coarseI, const IContainer& coarseIndexI, const int length)
  {
    const Array& coarse = dynamic_cast<const Array&>(coarseI);
    const Array<int>& coarseIndex = dynamic_cast<const Array<int>& >(coarseIndexI);

    for(int i=0; i<length; i++)
      _data[i] += coarse[coarseIndex[i]];
  }

  virtual shared_ptr<IContainer>
  newClone() const
  {
    return shared_ptr<Array>(new Array(_length));
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

#endif
