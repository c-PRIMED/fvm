// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _PRISM_H_
#define _PRISM_H_

#include "Cell.h"

struct Prism
{
public:

  enum{numFaces = 5};
  enum{numNodes = 6};

  static unsigned int faceNodeCount[numFaces];
  static unsigned int faceNodes[numFaces][4];
  
};


#endif
