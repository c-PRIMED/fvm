// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CELL_H_
#define _CELL_H_

#include <map>
using namespace std;
#include "Array.h"
#include "Vector.h"

class CRConnectivity;

template<class CellTrait>
class Cell
{
public:

  enum{numFaces = CellTrait::numFaces};
  enum{numNodes = CellTrait::numNodes};
  
  Cell();

  void
  orderCellFacesAndNodes(const int c,
                         CRConnectivity& cellFaces,
                         CRConnectivity& cellNodes,
                         const CRConnectivity& faceNodes,
                         const CRConnectivity& faceCells,
                         const Array<Vector<double,3> >& nodeCoordinates);

private:

  unsigned int _faceAllNodesSig[numFaces];
  unsigned int _nodeFirstFaceSig[numNodes];
  unsigned int _faceFirstFaceNodesSig[numNodes];

  map<unsigned int, unsigned int> _faceFirstFaceNodesSigMap;

};

// function to order connectivities over an entire mesh
void
orderCellFacesAndNodes(CRConnectivity& cellFaces,
                       CRConnectivity& cellNodes,
                       const CRConnectivity& faceNodes,
                       const CRConnectivity& faceCells,
                       const Array<Vector<double,3> >& nodeCoordinates);

#endif
