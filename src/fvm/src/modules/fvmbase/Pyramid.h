// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PYRAMID_H_
#define _PYRAMID_H_

#include "Cell.h"

struct Pyramid
{
public:

  enum{numFaces = 5};
  enum{numNodes = 5};

  static unsigned int faceNodeCount[numFaces];
  static unsigned int faceNodes[numFaces][4];
  
};


#endif
