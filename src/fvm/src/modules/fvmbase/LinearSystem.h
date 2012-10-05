// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _LINEARSYSTEM_H_
#define _LINEARSYSTEM_H_

#include "MultiFieldMatrix.h"


class LinearSystem
{
public:
 
  friend class LinearSystemMerger;

  LinearSystem();

  virtual ~LinearSystem();

  void initAssembly();

  void initSolve();

  void postSolve();

  void updateSolution();

  shared_ptr<LinearSystem>
  createCoarse(const int groupSize, const double weightRatioThreshold);

  MultiField& getX() {return *_x;}
  MultiField& getB() {return *_b;}
  MultiField& getDelta() {return *_delta;}
  MultiField& getResidual() {return *_residual;}
  
  MultiFieldMatrix& getMatrix() {return _matrix;}

  MultiField& getCoarseIndex() {return _coarseIndex;}

  shared_ptr<MultiField> getDeltaPtr() {return _delta;}
  shared_ptr<MultiField> getBPtr() {return _b;}

  void replaceDelta(shared_ptr<MultiField> newDelta) {_delta=newDelta;}
  void replaceB(shared_ptr<MultiField> newB) {_b=newB;}
  void replaceResidual(shared_ptr<MultiField> newR) {_residual=newR;}
  
  void setCoarseningField(const Field& f) {_coarseningField=&f;}

  bool isSymmetric;
  
private:
  MultiFieldMatrix _matrix;
  shared_ptr<MultiField> _x;
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
