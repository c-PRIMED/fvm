#ifndef _MATRIXJML_H_
#define _MATRIXJML_H_

template<class T>
class MatrixJML
{
 public:

  MatrixJML() {}
  virtual ~MatrixJML() {}
  virtual T& getElement(const int i,const int j)=0;
  virtual void zero()=0;

};

#endif
