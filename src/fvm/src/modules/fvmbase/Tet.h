#ifndef _TET_H_
#define _TET_H_

#include "Cell.h"

struct Tet
{
public:

  enum{numFaces = 4};
  enum{numNodes = 4};

  static unsigned int faceNodeCount[numFaces];
  static unsigned int faceNodes[numFaces][3];
  
};


#endif
