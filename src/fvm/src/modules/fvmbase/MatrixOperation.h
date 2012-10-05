// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MATRIXOPERATION_H_
#define _MATRIXOPERATION_H_

template <class T, int N>
class SquareMatrix
{
 public:
  enum { NSQR = N*N };

  SquareMatrix() {}
  ~SquareMatrix() {}

  SquareMatrix( const SquareMatrix& o){
    for ( int i=0; i<NSQR; i++){
      _data[i] = o._data[i];
    }
  }

  SquareMatrix(const T& s){
    for ( int i=0; i<N; i++){
      for ( int j=0; j<N; j++){
	if (i==j)
	  (*this)(i,j) = s;
	else
	  (*this)(i,j) = T(0);
      }
    }
  }

  T& operator()(int i, int j) {return _data[i*N+j];}

  SquareMatrix& operator=(const T& s)
  {
    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        (*this)(i,j) = (i==j) ? s : T(0.);
    return *this;
  }
  
  SquareMatrix& operator=(const SquareMatrix& o)
  {
    for(int i=0;i<NSQR;i++)
      _data[i] = o._data[i];
    return *this;
  }

  void print(){
    for(int i=0; i<N; i++){
      for (int j=0; j<N; j++)
	cout << (*this)(i,j) << " ";
      cout << endl;
    }
  }
 
 
 private:

  T _data[NSQR];

};

template<class T>
T determinant(SquareMatrix<T,2>& a){
  T d = a(0,0)*a(1,1)-a(0,1)*a(1,0);
  return d;
}

template<class T>
T determinant(SquareMatrix<T,3>& a){
  T d = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
    -a(0,1)*(a(1,0)*a(2,2) - a(1,2)*a(2,0))
    +a(0,2)*(a(1,0)*a(2,1) - a(1,1)*a(2,0));
  return d;
}

template<class T,int N>
T determinant(SquareMatrix<T,N>& a, int size){

  T d(0.0);
  T s(1.0);
    int m, n;
    SquareMatrix<T,N> b(NumTypeTraits<T>::getZero());
    
    if (size == 1){
      d = a(0,0);
      return d;
    }
    else{
      for(int k=0; k<size; k++){
	m = 0;
	n = 0;
	for(int i=0; i<size; i++){
	  for(int j=0; j<size; j++){
	    b(i,j)=0.0;
	    if(i!=0&&j!=k){
	      b(m,n)=a(i,j);
	      if(n<(size-2))
		n++;
	      else{
		n=0;
		m++;
	      }
	    }
	  }
	}
	d = d + s*(a(0,k)*determinant(b, size-1));
	s=-1*s;
      }      
    }  
    return d;
}

template<class T>
SquareMatrix<T,2> inverse(SquareMatrix<T,2>& a){

  SquareMatrix<T,2> inv;
  T d = determinant(a);
  inv(0,0) = a(1,1) / d;
  inv(0,1) = -a(0,1) / d;
  inv(1,0) = -a(1,0) / d;
  inv(1,1) = a(0,0) / d;
  
  return inv;
}


template<class T>
SquareMatrix<T,3>  inverse(SquareMatrix<T,3>& a)
{
  SquareMatrix<T,3> inv;
  T d = determinant(a);
  
  inv(0,0) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1)) / d;
  inv(0,1) =  (a(0,2)*a(2,1) - a(0,1)*a(2,2)) / d;
  inv(0,2) =  (a(0,1)*a(1,2) - a(0,2)*a(1,1)) / d;
  inv(1,0) =  (a(1,2)*a(2,0) - a(1,0)*a(2,2)) / d;
  inv(1,1) =  (a(0,0)*a(2,2) - a(0,2)*a(2,0)) / d;
  inv(1,2) =  (a(0,2)*a(1,0) - a(0,0)*a(1,2)) / d;
  inv(2,0) =  (a(1,0)*a(2,1) - a(1,1)*a(2,0)) / d;
  inv(2,1) =  (a(0,1)*a(2,0) - a(0,0)*a(2,1)) / d;
  inv(2,2) =  (a(0,0)*a(1,1) - a(0,1)*a(1,0)) / d;

  return inv;
}

template<class T, int N>
SquareMatrix<T,N>  inverse(SquareMatrix<T,N>& b, int size)
{
   const T d = determinant(b, size);
      
   SquareMatrix<T,N> a(NumTypeTraits<T>::getZero());
   SquareMatrix<T,N> tmp(NumTypeTraits<T>::getZero());
   SquareMatrix<T,N> fac(NumTypeTraits<T>::getZero());
   SquareMatrix<T,N> trans(NumTypeTraits<T>::getZero());
   
   int p, q, m, n, i, j;
   
   for(q=0; q<size; q++){
     for(p=0; p<size; p++){
       m=0;
       n=0;
       for(i=0; i<size; i++){
	 for(j=0; j<size; j++){
	   tmp(i,j)=0;
	   if(i!=q&&j!=p) {
	     tmp(m,n)=b(i,j);
	     if(n<(size-2))
	       n++;
	     else {
	       n=0;
	       m++;
	     }
	   }
	 }
       }
       fac(q,p)=pow(-1.0,q+p)*determinant(tmp,size-1);
     }
   }
  
   for(i=0; i<size; i++) {
     for(j=0; j<size; j++){
       trans(i,j)=fac(j,i);
     }
   }

   for(i=0; i<size; i++) {
     for(j=0; j<size; j++){
       a(i,j)=trans(i,j)/d;
     }
   }
   return a;
}
/* 
taken from http://www.physics.unlv.edu/~pang/cp_c.html 4.4
*/
template<class T, int N>
SquareMatrix<T,N> inverseGauss (SquareMatrix<T,N>&  a, int size)
{
  
  int indx[10];
  T c[10];
  SquareMatrix<T,N> b(NumTypeTraits<T>::getZero());
  SquareMatrix<T,N> x(NumTypeTraits<T>::getZero());
  
  b = 1.0; //identity matrix 
  /*
    for (int i = 0; i < size; i++)
    cout <<i<<" a-matrix "<<a(i,0)<< " "<< a(i,1)<< "   " <<a(i,2)<< " "<< a(i,3)<< "   " <<a(i,4)<<endl;		       
  */

 /* Initialize the index */
  for (int i = 0; i < size; i++){indx[i] = i;}
  
  /* Find the rescaling factors, one from each row */ 
  for (int i = 0; i < size; i++)
  {
    T c1 = 0;
    for (int j = 0; j < size; j++)
      {
	if (fabs(a(i,j)) > c1) 
	  c1 = fabs(a(i,j));
      }
    c[i] = c1;
  }
 
  
  /* Search the pivoting (largest) element from each column */ 
  for (int j = 0; j < size-1; j++)
  {
    T pi1 = 0.0; 
    int pk=j;
    int  itmp;
    T pi;
    for (int i = j; i < size; i++)
    {
      pi = fabs(a(indx[i],j))/c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        pk = i;
      }
    }
    /* Interchange the rows via indx[] to record pivoting order */
    itmp = indx[j];
    indx[j] = indx[pk];
    indx[pk] = itmp;
    for (int i = j+1; i < size; i++)
    {
      T pj;
      pj = a(indx[i],j)/a(indx[j],j);
      /* Record pivoting ratios below the diagonal */
      a(indx[i],j) = pj;
      /* Modify other elements accordingly */
      for (int k = j+1; k < size; k++)
      {
        a(indx[i],k) = a(indx[i],k)-pj*a(indx[j],k);
      }
    }
  }

  /*
  cout <<"indx  "<< "  re-scaling" <<endl;
  for (int i = 0; i < size; i++)cout<< indx[i]<<" -  "<<c[i] <<endl;
  for (int i = 0; i < size; i++)
    cout <<i<<" an-matrix "<<a(i,0)<< " "<< a(i,1)<< "   " <<a(i,2)<< " "<< a(i,3)<< "   " <<a(i,4)<<endl;	
  */

  for (int i = 0; i < size-1; i++)
  {
    for (int j = i+1; j < size; j++)
    {
      for (int k = 0; k < size; k++)
      {
        b(indx[j],k) = b(indx[j],k)-a(indx[j],i)*b(indx[i],k);
      }
    }
  }
  //backsubstitution
  for (int col = 0; col < size; col++)
  {
    x(size-1,col) = b(indx[size-1],col)/a(indx[size-1],size-1);
    for (int j = size-2; j >= 0; j = j-1)
    {
      x(j,col) = b(indx[j],col);
      for (int k = j+1; k < size; k++)
      {
        x(j,col) = x(j,col)-a(indx[j],k)*x(k,col);
      }
      x(j,col) = x(j,col)/a(indx[j],j);
    }
  }
  return x;
}

#endif
