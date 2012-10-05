// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _TRI_H_
#define _TRI_H_

#include "Cell.h"

struct Tri
{
public:

  enum{numFaces = 3};
  enum{numNodes = 3};

  static unsigned int faceNodeCount[numFaces];
  static unsigned int faceNodes[numFaces][2];
  
};


#endif
