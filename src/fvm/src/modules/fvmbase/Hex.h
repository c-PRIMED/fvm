// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _HEX_H_
#define _HEX_H_

#include "Cell.h"

struct Hex
{
public:

  enum{numFaces = 6};
  enum{numNodes = 8};

  static unsigned int faceNodeCount[numFaces];
  static unsigned int faceNodes[numFaces][4];
  
};


#endif
