// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "AABB.h"
#include "CRConnectivity.h"

#include <stack>

AABB::AABB(const Mesh& mesh)
{
  _is2D = mesh.getDimension() == 2;
  
  const Array<Vector<double,3> >& meshCoords = mesh.getNodeCoordinates();
  
  
  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
  {
      const FaceGroup& fg = *fgPtr;
      const StorageSite& faces = fg.site;
      const CRConnectivity& faceNodes = mesh.getFaceNodes(faces);
      
      const int nFaces = faces.getCount();
      
      for(int f=0; f<nFaces; f++)
      {
          if (_is2D)
          {
              _triangles.push_back(new MyTriangle(f,meshCoords,
                                                  faceNodes(f,0),
                                                  faceNodes(f,1)));
          }
          else
          {
              _triangles.push_back(new MyTriangle(f,0,meshCoords,
                                                  faceNodes(f,0),
                                                  faceNodes(f,1),
                                                  faceNodes(f,2)));
              if (faceNodes.getCount(f) == 4)
              {
                  
                  _triangles.push_back(new MyTriangle(f,1,meshCoords,
                                                      faceNodes(f,2),
                                                      faceNodes(f,3),
                                                      faceNodes(f,0)));
              }
          }
      }
  }
  
  if (_is2D)
    _tree_2D = boost::shared_ptr<CGAL_Tree_2D>(new CGAL_Tree_2D(_triangles.begin(),
                                                                _triangles.end()));
  else
  {
      _tree = boost::shared_ptr<CGAL_Tree>(new CGAL_Tree(_triangles.begin(),
                                                         _triangles.end()));
      
      //  _bbox = _tree->bbox();
  }
}

bool
AABB::hasIntersectionWithSegment(AABB::Vec3D a, AABB::Vec3D b)
{
  if (_is2D)
  {
      // can't do this check sine we are using 3d segments
      return false;
  }
  else
  {
      K::Segment_3 query(K::Point_3(a[0], a[1], a[2]),
                         K::Point_3(b[0], b[1], b[2]));
      return _tree->do_intersect(query); 
  }
}


bool
AABB::hasIntersectionWithTriangle(AABB::Vec3D a, AABB::Vec3D b, AABB::Vec3D c)
{
  K::Triangle_3 query(K::Point_3(a[0], a[1], a[2]),
                      K::Point_3(b[0], b[1], b[2]),
                      K::Point_3(c[0], c[1], c[2])
                      );
  if (_is2D)
    return _tree_2D->do_intersect(query);
  else
    return _tree->do_intersect(query); 
}

int
AABB::meshIntersections(const Mesh& mesh)
{
  const Array<Vector<double,3> >& meshCoords = mesh.getNodeCoordinates();
  int nIntersections = 0;

  if (_is2D)
  {
      const StorageSite& cells = mesh.getCells();
      const CRConnectivity& cellNodes = mesh.getCellNodes();
      
      const int nCells = cells.getSelfCount();
      for(int n=0; n<nCells; n++)
      {
          const Vec3D& a = meshCoords[cellNodes(n,0)];
          const Vec3D& b = meshCoords[cellNodes(n,1)];
          const Vec3D& c = meshCoords[cellNodes(n,2)];
          
          if (hasIntersectionWithTriangle(a,b,c))
          {
              nIntersections++;
          }
          else if (cellNodes.getCount(n) == 4)
          {
              
              const Vec3D& d = meshCoords[cellNodes(n,3)];
              if (hasIntersectionWithTriangle(c,d,a))
              {
                  nIntersections++;
              }
          }
      }
  }
  else
  {
      const StorageSite& faces = mesh.getFaces();
      const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
      
      const int nFaces = faces.getCount();
      for(int n=0; n<nFaces; n++)
      {
          const Vec3D& a = meshCoords[faceNodes(n,0)];
          const Vec3D& b = meshCoords[faceNodes(n,1)];
          const Vec3D& c = meshCoords[faceNodes(n,2)];
          
          if (hasIntersectionWithTriangle(a,b,c))
          {
              nIntersections++;
          }
          else if (faceNodes.getCount(n) == 4)
          {
              
              const Vec3D& d = meshCoords[faceNodes(n,3)];
              if (hasIntersectionWithTriangle(c,d,a))
              {
                  nIntersections++;
              }
          }
      }
  }
  return nIntersections;
}

int AABB::findOrientedSide(AABB::Vec3D p)
{
  if (_is2D)
  {
      K::Point_2 query(p[0], p[1]);
      foreach(const MyTriangle* myt, _triangles)
      {
          Line2D line = myt->getLine2D();
          CGAL::Oriented_side orientation = line.oriented_side(query);
          if (orientation == CGAL::ON_POSITIVE_SIDE)
            return 1;
          else if (orientation == CGAL::ON_ORIENTED_BOUNDARY)
          {
              K::Segment_2 segment = myt->getSegment2D();
              if (segment.has_on(query))
                return 0;
          }
      }
      return -1;
  }
  else
  {
      K::Point_3 query(p[0], p[1], p[2]);
      foreach(const MyTriangle* myt, _triangles)
      {
          Plane plane = myt->getPlane();
          CGAL::Oriented_side orientation = plane.oriented_side(query);
          if (orientation == CGAL::ON_POSITIVE_SIDE)
            return 1;
          else if (orientation == CGAL::ON_ORIENTED_BOUNDARY)
          {
              Triangle triangle = myt->getTriangle();
              if (triangle.has_on(query))
                return 0;
          }
      }
      return -1;
  }
}

