
#include "Mesh.h"
//#include "Field.h"
#include "StorageSite.h"

FieldLabel Mesh::coordinate("coordinate");
FieldLabel Mesh::area("area");
FieldLabel Mesh::areaMag("areaMag");
FieldLabel Mesh::volume("volume");

Mesh::Mesh(const int dimension):
  _dimension(dimension),
  _cells(0),
  _faces(0),
  _nodes(0),
  _interiorFaceGroup(),
  _faceGroups(),
  _boundaryGroups(),
  _interfaceGroups(),
  _connectivityMap()
{
  logCtor();
}

Mesh::~Mesh()
{
  logDtor();
}

const StorageSite&
Mesh::createInteriorFaceGroup(const int size)
{
  _interiorFaceGroup = shared_ptr<FaceGroup>(new FaceGroup(size,0,_faces,0,"interior"));
  _faceGroups.push_back(_interiorFaceGroup);
  return _interiorFaceGroup->site;
}


const StorageSite&
Mesh::createInterfaceGroup(const int size, const int offset, const int id)
{
  shared_ptr<FaceGroup> fg(new FaceGroup(size,offset,_faces,id,"interface"));
  _faceGroups.push_back(fg);
  _interfaceGroups.push_back(fg);
  return fg->site;
}


const StorageSite&
Mesh::createBoundaryFaceGroup(const int size, const int offset, const int id, const string& boundaryType)
{
  shared_ptr<FaceGroup> fg(new FaceGroup(size,offset,_faces,id,boundaryType));
  _faceGroups.push_back(fg);
  _boundaryGroups.push_back(fg);
  return fg->site;
}
