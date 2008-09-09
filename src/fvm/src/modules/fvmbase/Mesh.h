#ifndef _MESH_H_
#define _MESH_H_

#include "Array.h"

#include "StorageSite.h"
class CRConnectivity;
//#include "Field.h"
#include "FieldLabel.h"


class Mesh 
{
public:

  static FieldLabel coordinate;
  static FieldLabel area;
  static FieldLabel areaMag;
  static FieldLabel volume;
  
  typedef pair<const StorageSite*, const StorageSite*> SSPair;
  typedef map<SSPair,shared_ptr<CRConnectivity> > ConnectivityMap;
  enum
    {
      CELL_BAR2,
      CELL_TRI3,
      CELL_QUAD4,
      CELL_TETRA4,
      CELL_HEX8,
      CELL_PYRAMID,
      CELL_PRISM,
      CELL_TYPE_MAX
    } CellType;
  
  struct FaceGroup
  {
    FaceGroup(const int count_,
              const int offset_,
              const StorageSite& parent_,
              const int id_,
              const string& groupType_) :
      site(count_,0,offset_,&parent_),
      id(id_),
      groupType(groupType_)
    {}
    
    const StorageSite site;
    const int id;
    const string groupType;
  };

  typedef vector<shared_ptr<FaceGroup> > FaceGroupList;  

  Mesh(const int dimension);
  
  ~Mesh();

  DEFINE_TYPENAME("Mesh");

  const StorageSite& getFaces() const {return _faces;}
  const StorageSite& getCells() const {return _cells;}
  const StorageSite& getNodes() const {return _nodes;}

  StorageSite& getFaces() {return _faces;}
  StorageSite& getCells() {return _cells;}
  StorageSite& getNodes() {return _nodes;}

  //  const CRConnectivity& getConnectivity(const StorageSite& from,
  //                                      const StorageSite& to) const;

  //const CRConnectivity& getAllFaceCells() const;
  //const CRConnectivity& getFaceCells(const StorageSite& site) const;

  //const Array<int>& getCellTypes() const;
  //const Array<int>& getCellTypeCount() const;

  int getFaceGroupCount() const {return _faceGroups.size();}
  int getBoundaryGroupCount() const {return _boundaryGroups.size();}
  int getInterfaceGroupCount() const {return _interfaceGroups.size();}

  const FaceGroup&
  getFaceGroup(const int i) const {return *_faceGroups[i];}

  const FaceGroup&
  getBoundaryGroup(const int i) const {return *_boundaryGroups[i];}

  const FaceGroup&
  getInterfaceGroup(const int i) const {return *_interfaceGroups[i];}


  const StorageSite& createInteriorFaceGroup(const int size);
  const StorageSite& createInterfaceGroup(const int size,const int offset, 
                                    const int id);
  const StorageSite& createBoundaryFaceGroup(const int size, const int offset, 
                                       const int id, const string& boundaryType);
  
protected:
  const int _dimension;
  StorageSite _cells;
  StorageSite _faces;
  StorageSite _nodes;
  shared_ptr<FaceGroup> _interiorFaceGroup;
  FaceGroupList _faceGroups;
  FaceGroupList _boundaryGroups;
  FaceGroupList _interfaceGroups;
  ConnectivityMap _connectivityMap;
  map<const StorageSite*,const CRConnectivity*> _faceCellsMap;
  //mutable Array<int> *_cellTypes;
  //mutable Array<int> *_cellTypeCount;  
};

typedef vector<Mesh*> MeshList;

#endif
