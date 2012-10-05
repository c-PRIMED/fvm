// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _QUAD_H_
#define _QUAD_H_

#include "Cell.h"

struct Quad
{
public:

  enum{numFaces = 4};
  enum{numNodes = 4};

  static unsigned int getFaceNodeCount(const int) {return 2;}

  static unsigned int faceNodeCount[numFaces];
  static unsigned int faceNodes[numFaces][2];
  
};


#endif
