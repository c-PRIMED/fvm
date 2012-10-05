// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _ARRAY2D_H_
#define _ARRAY2D_H_

#include "ArrayBase.h"
#include <iostream>
#include <iomanip>
#include "NumType.h"

template <class T>
class Array2D
{

public:

  
  explicit Array2D(const int row_size, const int col_size):
    _self(*this),
    _rowSize(row_size),
    _colSize(col_size),
    _length(row_size*col_size),
    _data(new T[_length])
  {
    init();
  }
  
  
   ~Array2D()
  {
     delete [] _data;
  }

  
  int getRow() const { return _rowSize; }
  int getCol() const { return _colSize; }


  /**
   * returns the index of the j'th non zero column for row i
   *
   */

  T& operator()(const int i, const int j)
  {
    return _data[i*_colSize+j];
  }

  const T& operator()(const int i, const int j) const
  {
    return _data[i*_colSize+j];
  }
  void operator=(const T& x){
     for ( int i = 0; i < _length; i++ )
         _data[i] = x;
  }        
  //copy of aij to this array
  void partialCopyFrom(const Array2D& aij){
      const int nrow = aij.getRow();
      const int ncol = aij.getCol();
      for ( int i = 0; i < nrow; i++ ){
          for ( int j = 0; j < ncol; j++){
	       _self(i,j) = aij(i,j); 
	  }
      }
  }
  //partial copy of this array to aij 
  void partialCopyTo(Array2D& aij){
      const int nrow = aij.getRow();
      const int ncol = aij.getCol();
      for ( int i = 0; i < nrow; i++ ){
          for ( int j = 0; j < ncol; j++){
	       aij(i,j) = _self(i,j);   
	  }
      }
  }
 
  void zeros()
  {
     for ( int i = 0; i < _length; i++ )
         _data[i] =  NumTypeTraits<T>::getZero();
  }
  void setIdentity()
  {
      for ( int i = 0; i < _rowSize; i++ ){
          for ( int j = 0; j < _colSize; j++){
	      if ( i == j )
	         _self(i,j) = NumTypeTraits<T>::getUnity();

	  }
      }
  }


  void*  getData() { return _data;} 
  int getDataSize() const
  {
    return  _length*NumTypeTraits<T>::getDataSize();
  }


 void print(ostream& os) const
  {
      for ( int i = 0; i < _rowSize; i++){
         for ( int j = 0; j < _colSize; j++ ){
            os << std::setprecision(14) << this->operator()(i,j) << "    ";
	 }
	 os << "\n";
      }
      os << endl;
  }

private:
  Array2D(const Array2D&);
  void init()
  {
    for ( int i = 0; i < _length; i++ )
        _data[i] =  NumTypeTraits<T>::getZero();//NumTypeTraits<T>::getNegativeUnity();
  }

  Array2D&   _self;

  int _rowSize;
  int _colSize;
  int _length;
  T* _data;
};


#endif
