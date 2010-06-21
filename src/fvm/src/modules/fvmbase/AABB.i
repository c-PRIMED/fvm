%{
#include "AABB.h"
  %}

class AABB
{
public:
  typedef Vector<double,3> Vec3D;
  AABB(const Mesh& mesh);
  bool hasIntersectionWithSegment(Vec3D a, Vec3D b);
  bool hasIntersectionWithTriangle(Vec3D a, Vec3D b, Vec3D c);
  int meshIntersections(const Mesh& mesh);
  int findOrientedSide(Vec3D p);
};
