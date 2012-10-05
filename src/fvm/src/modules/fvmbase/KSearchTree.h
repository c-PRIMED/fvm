// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KSEARCHTREE_H_
#define _KSEARCHTREE_H_

#include "GeomFields.h"

#include "Mesh.h"

#include <CGAL/basic.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>



struct MyPoint
{
  typedef Vector<double,3> Vec3D;

  Vec3D vec;
  int index;

  MyPoint(const MyPoint& o) :
    vec(o.vec),
    index(o.index)
  {}
  
  MyPoint(const Vec3D& vec_,  const int index_) :
    vec(vec_),
    index(index_)
  {}
};

namespace CGAL
{

template <>
struct Kernel_traits<MyPoint> {
  struct Kernel {
    typedef double FT;
    typedef double RT;
  };
};
}

struct Construct_coord_iterator
{
  const double* operator()(const MyPoint& p) const
  { return &(p.vec[0]); }

  const double* operator()(const MyPoint& p, int)  const
  { return &(p.vec[0])+3; }
};


struct Distance
{
  typedef MyPoint Query_item;
  
  double transformed_distance(const MyPoint& p1, const MyPoint& p2) const
  {
    double distx= p1.vec[0]-p2.vec[0];
    double disty= p1.vec[1]-p2.vec[1];
    double distz= p1.vec[2]-p2.vec[2];
    return distx*distx+disty*disty+distz*distz;
  }
  
  template <class TreeTraits>
  double min_distance_to_rectangle(const MyPoint& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const
  {
    double distance(0.0), h = p.vec[0];
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.vec[1];
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.vec[2];
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }

  template <class TreeTraits>
  double max_distance_to_rectangle(const MyPoint& p,
				   const CGAL::Kd_tree_rectangle<TreeTraits>& b) const
  {
    double h = p.vec[0]; 	  
    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
      (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);
    
    h=p.vec[1];
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
      (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);

    h=p.vec[2];
    double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
      (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }
  
  double new_distance(double& dist, double old_off, double new_off,
		      int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }
  
  double transformed_distance(double d) const { return d*d; }
  
  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }
  
}; // end of struct Distance


/**
 * A wrapper for CGAL's K_neighbor_search
 * 
 */

class KSearchTree
{
public:
  typedef Vector<double,3> Vec3D;
  typedef Array<Vec3D> Vec3DArray;
  typedef Array<int> IntArray;

  KSearchTree();
  KSearchTree(const Vec3DArray& points);
  void insert(const Vec3D& v, const int n);

  void findNeighbors(const Vec3D& p, const int k, Array<int>& neighbors);
  
private:


  typedef CGAL::Search_traits<double, MyPoint, const double*,
                              Construct_coord_iterator> Traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Traits> K_neighbor_search;
  typedef K_neighbor_search::Tree Tree;

  boost::shared_ptr<Tree> _tree;
};

#endif
