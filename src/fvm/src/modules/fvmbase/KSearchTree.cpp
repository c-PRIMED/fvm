// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "KSearchTree.h"


KSearchTree::KSearchTree(const Vec3DArray& points)
{
  _tree = boost::shared_ptr<Tree>(new Tree());

  const int nPoints = points.getLength();
  for(int n=0; n<nPoints; n++)
    _tree->insert(MyPoint(points[n],n));
}

KSearchTree::KSearchTree()
{
  _tree = boost::shared_ptr<Tree>(new Tree());
}

void
KSearchTree::insert(const Vec3D& v, const int n)
{
  _tree->insert(MyPoint(v,n));
}


void
KSearchTree::findNeighbors(const Vec3D& p, const int k, Array<int>& neighbors)
{
  MyPoint query(p, -1);
  if ( _tree->size() == 0 ){
     return;
  }
  K_neighbor_search search(*_tree, query, k);

  int i=0;
  for(K_neighbor_search::iterator it = search.begin();
      it != search.end();
      it++)
  {
      neighbors[i++] = it->first.index;
  }

}
