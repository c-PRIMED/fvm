// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _AABB_H_
#define _AABB_H_

#include "GeomFields.h"

#include "Mesh.h"

#include <CGAL/AABB_tree.h> // must be inserted before kernel
#include <CGAL/AABB_traits.h>

#include <CGAL/Simple_cartesian.h>



/**
 * A wrapper for CGAL's axis aligned bounding box concept. We use this
 * to compute intersections of the fluid mesh with the boundary mesh
 * representing the solid. The tree is constructed with the faces of
 * the latter and then we go through all the faces of the fluid mesh
 * to find the ones that have intersections with any of the tree
 * faces.
 * 
 */

class AABB
{
public:
  typedef Vector<double,3> Vec3D;

  AABB(const Mesh& mesh);

  /**
   * check whether a segment defined by the two points a and b has any
   * intersections with the faces of the tree
   * 
   */

  bool hasIntersectionWithSegment(Vec3D a, Vec3D b);

    /**
   * check whether a triangle defined by the three points a, b and c has any
   * intersections with the faces of the tree
   * 
   */

  bool hasIntersectionWithTriangle(Vec3D a, Vec3D b, Vec3D c);

  
  int meshIntersections(const Mesh& mesh);

  /**
   * locate the side on which a given point lies. returns 0 if the
   * point is on one of the faces of the tree, -1 if it's inside and 1
   * if it is outside. This check is only meaningful if the faces of
   * the tree form a closed surface.
   * 
   */

  int findOrientedSide(Vec3D p);
  
private:

  typedef CGAL::Simple_cartesian<double> K;
  typedef K::Point_3    Point; // CGAL 3D point type
  typedef K::Point_2    Point2D; // CGAL 2D point type

  typedef K::Triangle_3    Triangle; // CGAL 3D triangle type
  typedef K::Plane_3    Plane; // CGAL 3D plane
  typedef K::Line_2    Line2D; // CGAL 2d line

  //typedef CGAL::Bbox_3    BoundingBox;

  
  struct MyTriangle
  {
    const int faceIndex;
    const int subFaceIndex;
    const Array<Vec3D>& coordArray;
    int vertices[3];
    
    MyTriangle(const int faceIndex_,
               const int subFaceIndex_,
               const Array<Vec3D>& coordArray_,
               const int v0_,
               const int v1_,
               const int v2_) :
      faceIndex(faceIndex_),
      subFaceIndex(subFaceIndex_),
      coordArray(coordArray_)
    {
      vertices[0] = v0_;
      vertices[1] = v1_;
      vertices[2] = v2_;
    }

    MyTriangle(const int faceIndex_,
               const Array<Vec3D>& coordArray_,
               const int v0_,
               const int v1_) :
      faceIndex(faceIndex_),
      subFaceIndex(0),
      coordArray(coordArray_)
    {
      vertices[0] = v0_;
      vertices[1] = v1_;
      vertices[2] = -1;
    }

    Point getVertex(const int n) const
    {
      const Vec3D& v = coordArray[vertices[n]];
      return Point(v[0], v[1], v[2]);
    }

    Point2D getVertex2D(const int n) const
    {
      const Vec3D& v = coordArray[vertices[n]];
      return Point2D(v[0], v[1]);
    }

    Triangle getTriangle() const
    {
      return Triangle(getVertex(0),
                      getVertex(1),
                      getVertex(2));
    }

    K::Segment_2 getSegment2D() const
    {
      return K::Segment_2(getVertex2D(0),
                          getVertex2D(1));
    }

    K::Segment_3 getSegment() const
    {
      return K::Segment_3(getVertex(0),
                          getVertex(1));
    }

    Plane getPlane() const
    {
      return Plane(getVertex(0),
                   getVertex(1),
                   getVertex(2));
    }
    
    Line2D getLine2D() const
    {
      // CGAL's line orientation is opposite ours in 2d
      return Line2D(getVertex2D(1),
                   getVertex2D(0));
    }
    
  };
  
  // the custom triangles are stored into a vector
  typedef std::vector<MyTriangle*>::const_iterator MyTriangleIterator;
  
  // The following primitive provides the conversion facilities between
  // the custom triangle and point types and the CGAL ones
  struct MyTrianglePrimitive
  {
  public:
    
    // this is the type of data that the queries returns
    
    typedef const MyTriangle* Id;
    
    // CGAL types returned
    typedef K::Triangle_3 Datum; // CGAL 3D triangle type
    
    
    // default constructor needed
    MyTrianglePrimitive()
    {} 

    MyTrianglePrimitive(MyTriangle* t) :
      m_pt(t)
    {} 
 
    // the following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree
    MyTrianglePrimitive(MyTriangleIterator it)
      : m_pt(*it)
    {}
    
    const Id& id() const { return m_pt; }
  
    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
      return m_pt->getTriangle();
    }
    
    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return m_pt->getVertex(0);}
    
  private:
    Id m_pt; // this is what the AABB tree stores internally
    
  };
  
  // The following primitive provides the conversion facilities between
  // the custom triangle and point types and the CGAL ones
  struct MySegmentPrimitive
  {
  public:
    
    // this is the type of data that the queries returns
    
    typedef const MyTriangle* Id;
    
    // CGAL types returned
    typedef K::Segment_3 Datum; 
    
    
    // default constructor needed
    MySegmentPrimitive()
    {} 

    MySegmentPrimitive(MyTriangle* t) :
      m_pt(t)
    {} 
 
    // the following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree
    MySegmentPrimitive(MyTriangleIterator it)
      : m_pt(*it)
    {}
    
    const Id& id() const { return m_pt; }
  
    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
      return m_pt->getSegment();
    }
    
    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return m_pt->getVertex(0);}
    
  private:
    Id m_pt; // this is what the AABB tree stores internally
    
  };
  

  typedef CGAL::AABB_traits<K, MyTrianglePrimitive> My_AABB_traits;
  typedef CGAL::AABB_tree<My_AABB_traits> CGAL_Tree;
  
  typedef CGAL::AABB_traits<K, MySegmentPrimitive> My_AABB_traits_2D;
  typedef CGAL::AABB_tree<My_AABB_traits_2D> CGAL_Tree_2D;
  
  bool _is2D;
  std::vector<MyTriangle*> _triangles;
  boost::shared_ptr<CGAL_Tree> _tree;
  boost::shared_ptr<CGAL_Tree_2D> _tree_2D;
  //BoundingBox _bbox;
};

#endif
