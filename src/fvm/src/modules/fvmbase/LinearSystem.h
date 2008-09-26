#ifndef _LINEARSYSTEM_H_
#define _LINEARSYSTEM_H_

#include "MultiFieldMatrix.h"

class LinearSystem
{
public:
  LinearSystem();

  virtual ~LinearSystem();

  void initAssembly();

  void initSolve();

  void postSolve();

  void updateSolution();

  shared_ptr<LinearSystem>
  createCoarse(const int groupSize, const double weightRatioThreshold);

  MultiField& getX() {return _x;}
  MultiField& getB() {return *_b;}
  MultiField& getDelta() {return *_delta;}
  MultiField& getResidual() {return *_residual;}
  
  MultiFieldMatrix& getMatrix() {return _matrix;}

  MultiField& getCoarseIndex() {return _coarseIndex;}
  
private:
  MultiFieldMatrix _matrix;
  MultiField _x;
  shared_ptr<MultiField> _b;
  shared_ptr<MultiField> _delta;
  shared_ptr<MultiField> _residual;
  MultiField _coarseIndex;
  const Field* _coarseningField;
  shared_ptr<MultiField> _xAux;
  shared_ptr<MultiField> _bAux;
  shared_ptr<MultiField> _deltaAux;
  shared_ptr<MultiField> _residualAux;
};
#endif
