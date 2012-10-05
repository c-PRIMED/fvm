// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Cell.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include "CRConnectivity.h"

using namespace std;

#include "Quad.h"
#include "Tri.h"
#include "Hex.h"
#include "Tet.h"
#include "Pyramid.h"
#include "Prism.h"

unsigned int Quad::faceNodeCount[4] = {2,2,2,2};
unsigned int Quad::faceNodes[4][2] = { {0,1}, {1,2}, {2,3}, {3, 0}};

unsigned int Tri::faceNodeCount[3] = {2,2,2};
unsigned int Tri::faceNodes[3][2] = { {0,1}, {1,2}, {2,0}};

unsigned int Hex::faceNodeCount[6] = {4,4,4,4,4,4};
unsigned int Hex::faceNodes[6][4] = { {0,1,2,3},
                                      {4,7,6,5},
                                      {0,4,5,1},
                                      {1,5,6,2},
                                      {2,6,7,3},
                                      {3,7,4,0}};

unsigned int Tet::faceNodeCount[4] = {3,3,3,3};
unsigned int Tet::faceNodes[4][3] = { {0,1,2}, {0,3,1}, {1,3,2}, {2,3,0}};

unsigned int Pyramid::faceNodeCount[5] = {4,3,3,3,3};
unsigned int Pyramid::faceNodes[5][4] = { {0,1,2,3},
                                          {0,4,1,0},
                                          {1,4,2,0},
                                          {2,4,3,0},
                                          {3,4,0,0}};

unsigned int Prism::faceNodeCount[5] = {3,3,4,4,4};
unsigned int Prism::faceNodes[5][4] = { {0,1,2,0},
                                        {3,5,4,0},
                                        {0,3,4,1},
                                        {1,4,5,2},
                                        {2,5,3,0}};

template<class CellTrait>
Cell<CellTrait>::Cell()
{
      for(int n=0; n<numNodes; n++) _nodeFirstFaceSig[n]  = 0;

    for(int f=0; f<numFaces; f++)
    {
        _faceAllNodesSig[f] = 0;
        for(unsigned int nf=0; nf<CellTrait::faceNodeCount[f]; nf++)
        {
            unsigned int n = CellTrait::faceNodes[f][nf];
            _faceAllNodesSig[f] |= (1<<n);
            if (f == 0)
              _nodeFirstFaceSig[n] = 1<<n;
        }
    }

    for(int f=0; f<numFaces; f++)
    {
        _faceFirstFaceNodesSig[f] = 0;
        for(unsigned int nf=0; nf<CellTrait::faceNodeCount[f]; nf++)
        {
            unsigned int n = CellTrait::faceNodes[f][nf];
            _faceFirstFaceNodesSig[f] |= _nodeFirstFaceSig[n];
        }
      _faceFirstFaceNodesSigMap[_faceFirstFaceNodesSig[f]] = f;
    }

#if 0
    for(int f=0; f<numFaces; f++)
      cout << _faceAllNodesSig[f] << " ";
    cout << endl;
    
    for(int n=0; n<numNodes; n++)
      cout << _nodeFirstFaceSig[n] <<  " " ;
    cout << endl;

    for(int f=0; f<numFaces; f++)
      cout << _faceFirstFaceNodesSig[f] <<  "   " ;
    cout << endl;
#endif
}


template<class CellTrait>
void
Cell<CellTrait>::orderCellFacesAndNodes(const int c,
                                        CRConnectivity& cellFaces,
                                        CRConnectivity& cellNodes,
                                        const CRConnectivity& faceNodes,
                                        const CRConnectivity& faceCells,
                                        const Array<Vector<double,3> >& nodeCoordinates)
{

  int face0 = 0;
  int face0NodeCount=0;
  bool reverseFace0Nodes = false;

  for(int nf=0; nf<numFaces; nf++)
  {
      const int f = cellFaces(c,nf);
      if ((unsigned int)faceNodes.getCount(f) == CellTrait::faceNodeCount[0])
      {
          face0 = f;
          face0NodeCount = faceNodes.getCount(f);
          if (faceCells(f,0) != c)
          {
              if (faceCells(f,1) != c)
                cerr  << "malformed grid " << endl;
              reverseFace0Nodes = true;
          }
          break;
      }
  }

  map<int, unsigned int> thisNodesFirstFaceNodesSig;

  for(unsigned int nn=0; nn<(unsigned int)face0NodeCount; nn++)
  {
      int node = faceNodes(face0,nn);
      if (reverseFace0Nodes) node = faceNodes(face0,face0NodeCount-nn-1);
      thisNodesFirstFaceNodesSig[node] = 1<<nn;
  }

  int orderedFaces[numFaces];

  for(int nf=0; nf<numFaces; nf++)
  {
      unsigned int thisFaceFirstFaceNodesSig = 0;
      const int f = cellFaces(c,nf);

      for(int nn=0; nn<faceNodes.getCount(f); nn++)
      {
          const int n = faceNodes(f,nn);

          if (thisNodesFirstFaceNodesSig.find(n) !=
              thisNodesFirstFaceNodesSig.end())
            thisFaceFirstFaceNodesSig |= thisNodesFirstFaceNodesSig[n];
      }

      int faceOrder = _faceFirstFaceNodesSigMap[thisFaceFirstFaceNodesSig];
      orderedFaces[faceOrder] = f;
  }

  for(int nf=0; nf<numFaces; nf++)
    cellFaces(c,nf) = orderedFaces[nf];
  

  map<int, unsigned int> thisNodesAllFaceNodesSig;

  for(int nf=0; nf<numFaces; nf++)
  {
      //      unsigned int thisFaceFirstFaceNodeSig = 0;
      const int f = cellFaces(c,nf);

      //      bool reverseFaceNodes = (faceCells(f,0) != c);

      const int faceNodeCount = faceNodes.getCount(f);
      for(int nn=0; nn<faceNodeCount; nn++)
      {
          int n = faceNodes(f,nn);
          if (reverseFace0Nodes) n = faceNodes(f,faceNodeCount-nn-1);

          if (thisNodesAllFaceNodesSig.find(n) ==
              thisNodesAllFaceNodesSig.end())
          {
              thisNodesAllFaceNodesSig[n] =
                _faceAllNodesSig[nf];
          }
          else
          {
              thisNodesAllFaceNodesSig[n] &=
                _faceAllNodesSig[nf];
          }
      }
  }

  for( map<int,unsigned int>::const_iterator pos =
         thisNodesAllFaceNodesSig.begin();
       pos != thisNodesAllFaceNodesSig.end();
       ++pos)
  {
    const int node = pos->first;
    // log2 missing in msvc ?
          const int index = (int) log2(pos->second);
    // const int index = (int) (log( (double) pos->second)/log(2.));
      cellNodes(c,index) = node;
  }
   
}

template class Cell<Quad>;
template class Cell<Tri>;
template class Cell<Hex>;
template class Cell<Tet>;
template class Cell<Pyramid>;
template class Cell<Prism>;


static Cell<Quad> quad;
static Cell<Tri> tri;
static Cell<Hex> hexCell;
static Cell<Tet> tetCell;
static Cell<Pyramid> pyramidCell;
static Cell<Prism> prismCell;

#include "StorageSite.h"

struct MyCoords
{
  MyCoords(int i_, const double x_) :
    i(i_),
    x(x_)
  {}

  MyCoords(const MyCoords& o) :
    i(o.i),
    x(o.x)
  {}
  
  int i;
  double x;
};

bool myCoordComparison (const MyCoords& a, const MyCoords& b)
{
  return (a.x<b.x);
}

void
orderCellFacesAndNodes(CRConnectivity& cellFaces,
                       CRConnectivity& cellNodes,
                       const CRConnectivity& faceNodes,
                       const CRConnectivity& faceCells,
                       const Array<Vector<double,3> >& nodeCoordinates)
{

  const int numCells = cellNodes.getRowSite().getSelfCount();

  for(int c=0; c<numCells; c++)
  {
      const int numCellNodes = cellNodes.getCount(c);
      const int numCellFaces = cellFaces.getCount(c);

      int edgeFaceCount  =0;
      int triFaceCount = 0;
      int quadFaceCount = 0;
      for(int nf=0; nf<numCellFaces; nf++)
      {
          const int f = cellFaces(c,nf);
          const int faceNodeCount = faceNodes.getCount(f);
          switch(faceNodeCount)
          {
          case 2:
            edgeFaceCount++;
            break;
          case 3:
            triFaceCount++;
            break;
          case 4:
            quadFaceCount++;
            break;
          default:
            break;
          }
      }

      if (numCellNodes == 4 && edgeFaceCount == 4)
      {
          quad.orderCellFacesAndNodes(c,cellFaces,cellNodes,
                                      faceNodes,faceCells,nodeCoordinates);
      }
      else if (numCellNodes == 3 && edgeFaceCount == 3)
      {
          tri.orderCellFacesAndNodes(c,cellFaces,cellNodes,
                                     faceNodes,faceCells,nodeCoordinates);
      }
      else if (numCellNodes == edgeFaceCount)
      {
          vector<MyCoords> angles;
          Vector<double,3> mean(Vector<double,3>::getZero());
          for(int nn=0; nn<numCellNodes; nn++)
            mean += nodeCoordinates[cellNodes(c,nn)];
          mean /= numCellNodes;
          for(int nn=0; nn<numCellNodes; nn++)
          {
              const int nodeIndex = cellNodes(c,nn);
              const Vector<double,3> xm =
                nodeCoordinates[nodeIndex] - mean;
              const double angle = atan2(xm[1], xm[0]);
              angles.push_back(MyCoords(nodeIndex,angle));
          }

          sort(angles.begin(),angles.end(),myCoordComparison);
          
          for(int nn=0; nn<numCellNodes; nn++)
          {
              cellNodes(c,nn) = angles[nn].i;
          }
      }
      else if (numCellNodes == 8 && quadFaceCount == 6)
      {
          hexCell.orderCellFacesAndNodes(c,cellFaces,cellNodes,
                                         faceNodes,faceCells,nodeCoordinates);
      }
      else if (numCellNodes == 4 && triFaceCount == 4)
      {
          tetCell.orderCellFacesAndNodes(c,cellFaces,cellNodes,
                                         faceNodes,faceCells,nodeCoordinates);
      }
      else if (numCellNodes == 5 && triFaceCount == 4 &&
               quadFaceCount == 1)
      {
          pyramidCell.orderCellFacesAndNodes(c,cellFaces,cellNodes,
                                             faceNodes,faceCells,
                                             nodeCoordinates);
      }
      else if (numCellNodes == 6 && triFaceCount == 2 &&
               quadFaceCount == 3)
      {
          prismCell.orderCellFacesAndNodes(c,cellFaces,cellNodes,
                                           faceNodes,faceCells,
                                           nodeCoordinates);
      }
      else
      {
#if 0
          cerr << "unimplemented cell type: numCellnodes = "
               << numCellNodes
               << " quadFaceCount = " << quadFaceCount
               << " triFaceCount = " << triFaceCount
               << " edgeFaceCount = " << edgeFaceCount
               << endl;
#endif
      }
  }
}

