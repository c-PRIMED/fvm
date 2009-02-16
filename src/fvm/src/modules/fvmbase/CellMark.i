%{
#include "CellMark_impl.h"
#include "Field.h"
#include "Mesh.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "Vector.h"
#include "NumType.h"
#include "Octree.h"
#include "CellMark.h"
#include "Array.h"
#include "MPM_Particles.h"

  %}
%include "Vector.i"
%include "Mesh.i"
%include "Octree.i"
%include "MPM_Particles.i"

using namespace std;

void CellMark_Impl(Mesh& mesh, const GeomFields& geomFields, const string fileBase, Octree& O, MPM& solid);

 


