// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ARRAYWRITER_H_
#define _ARRAYWRITER_H_


#include "Array.h"
#include "Vector.h"

class ArrayWriter
{
public:
  ArrayWriter(const bool binary, const int vectorComponent,
              const int atypeComponent):
    _binary(binary),
    _vectorComponent(vectorComponent),
    _atypeComponent(atypeComponent)
  {}

  virtual ~ArrayWriter() {};
    
  template<class T>
  void writeFloats(FILE *fp, const T *data, const int count, const int stride,
                   const Array<bool> *mask=0) 
  {
    if (mask)
    {      
        for(int n=0,i=0; n<count; i+=stride,n++)
          if ((*mask)[n])
            fprintf(fp,"%12.5e\n",data[i]);
    }
    else
    {
        for(int n=0,i=0; n<count; i+=stride,n++)
          fprintf(fp,"%12.5e\n",data[i]);
    }
  }

protected:
  const bool _binary;
  const int _vectorComponent;
  const int _atypeComponent;
};
  
template<class T>
class ScalarArrayWriter : public ArrayWriter
{
public:
  typedef T ElementType;
  typedef Array<ElementType> ArrayType;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  
  ScalarArrayWriter(const bool binary, const int vectorComponent,
                    const int atypeComponent) :
    ArrayWriter(binary,vectorComponent,atypeComponent)
  {}
  
  void write(FILE *fp, const Array<T>& array,
             const int iBeg, int count,
             const Array<bool> *mask=0)
  {
    const int length = array.getLength();
    if (count == -1)
      count = length;
    else if (count+iBeg > length)
      throw CException("invalid count ");
    
    ElementType *dataArray = (ElementType*)(array.getData());
    
    T_BuiltIn *data = (T_BuiltIn*)(dataArray+iBeg);
    
    int stride = 1;
    
#ifdef USING_ATYPE_TANGENT
    if ((_atypeComponent == 0) || (_atypeComponent == 1))
    {
        data += _atypeComponent;
        stride = 2;
    }
    else
      throw CException("invalid component for Tangent");
        
#endif
    writeFloats(fp,data,count,stride,mask);
  }
};


template<class T, int N>
class VectorArrayWriter : public ArrayWriter
{
public:
  typedef Vector<T,N> ElementType;
  typedef Array<ElementType> ArrayType;
  typedef typename NumTypeTraits<T>::T_BuiltIn T_BuiltIn;
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  
  VectorArrayWriter(const bool binary, const int vectorComponent,
                    const int atypeComponent) :
    ArrayWriter(binary,vectorComponent,atypeComponent)
  {}

  virtual void write(FILE *fp, const ArrayType& array,
                     const int iBeg, int count,
                     const Array<bool> *mask=0)
  {
    const int length = array.getLength();
    if (count == -1)
      count = length;
    else if (count+iBeg > length)
      throw CException("invalid count ");
    
    
    ElementType *dataArray = (ElementType*)(array.getData());
    T_Scalar *dataScalar = (T_Scalar*)(dataArray+iBeg);
    
    int stride = 1;

    if (_vectorComponent == -2)
    {
        T_BuiltIn *data = (T_BuiltIn*)(dataScalar);

        // write all components in coupled order
        
#ifdef USING_ATYPE_TANGENT
        if ((_atypeComponent == 0) || (_atypeComponent == 1))
        {
            data += _atypeComponent;
            stride = 2;
        }
        else
          throw CException("invalid component for Tangent");
#endif
        if (mask)
          throw CException("masked write not implemented");
        writeFloats(fp,data,count*N,stride,0);
    }
    else 
    {
        // write all components in sequential order
        int ndBeg=0;
        int ndEnd=N;

        if (_vectorComponent >= 0 && _vectorComponent <N)
        {
            // restrict to just one component
            ndBeg=_vectorComponent;
            ndEnd=ndBeg+1;
        }
        else if (_vectorComponent != -1)
        {
            throw CException("invalid vector component");
        }
        
        stride = N;

        for(int nd=ndBeg; nd<ndEnd; nd++)
        {
            T_BuiltIn *nddata = (T_BuiltIn*) (dataScalar+nd);
            
#ifdef USING_ATYPE_TANGENT
            if ((_atypeComponent == 0) || (_atypeComponent == 1))
            {
                nddata += _atypeComponent;
                stride = N*2;
            }
            else
              throw CException("invalid component for Tangent");
#endif
            writeFloats(fp,nddata,count,stride,mask);
        }
    }   
  }
};


#endif
