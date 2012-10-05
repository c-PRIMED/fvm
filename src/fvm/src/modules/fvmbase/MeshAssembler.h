// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _MESHASSEMBLER_H_
#define _MESHASSEMBLER_H_

#include "Mesh.h"
#include "Array.h"
#include "Vector.h"
#include <vector>
#include <set>
#include <fstream>

class MeshAssembler
{
public:
    
    typedef   shared_ptr< Array<double> >  ArrayDblePtr;
    typedef   shared_ptr< Array<int> >     ArrayIntPtr;
    typedef   map<int,ArrayIntPtr>         ArrayIntPtrMap;
    typedef   shared_ptr< StorageSite >  StorageSitePtr; 
    typedef   shared_ptr< CRConnectivity > CRConnectivityPtr;
    //typedef   shared_ptr< Array<int> >     ArrayIntPtr;
    typedef   shared_ptr< Array<Mesh::VecD3> >   ArrayVecD3Ptr;
    typedef   vector< map<int, set<int> > >      VecMap;
    


   MeshAssembler( const MeshList& meshList );
   ~MeshAssembler();

  //return single mesh, for sake of compatibality we return as meshList
   const MeshList&  meshList() const { return _mesh;}

   void  debug_print();
   void  debug_sites();
   void  debug_localToGlobal_mappers();
   void  debug_globalCellToMeshID_mappers();
   void  debug_sync_localToGlobal_mappers();
   void  debug_faceCells();
   void  debug_localNodeToGlobal();

   //get methods


private:
   MeshAssembler(const MeshAssembler&);
   void  init();

   void setCellsSite();
   void setFacesSite();
   void setInterfaceNodes();
   void setNodesSite();
   void setCellsMapper(); 
   void setNodesMapper();
   void setBoundaryFaceGroup();
   void setCoord();
   void setSites();
   void setMeshCellColor();


   void setFaceCells();
   void setFaceNodes();
   void setMesh();  

   void countInterfaceNodes();
   
   int  getInnerNodesCount();
   int  getInterfaceNodesDuplicatedCount();
   int  getInterfaceNodesCount();

   void  debug_file_open( const string& fname );
   void  debug_file_close();

   const MeshList _meshList;

   StorageSitePtr _cellSite;
   StorageSitePtr _faceSite;
   StorageSitePtr _nodeSite;
  
   VecMap    _interfaceNodesSet; //[meshid][faceid] gives set for nodes
  //mappers for inner cells
   map<int, ArrayIntPtr >  _localCellToGlobal; //this include ghost cells at interfaces, ghostcells on boundary is set =-1
   vector<int>             _globalCellToMeshID;  //belongs to which mesh from global Cell ID
   vector<int>             _globalCellToLocal;   //gives local Cell ID from Global ID but which mesh comes from mapGlobalCellToMeshID

   vector < map<int,int> >  _localInterfaceNodesToGlobalMap; // for each mesh
   map<int, ArrayIntPtr >   _localNodeToGlobal; 

   set<Mesh::VecD3>  _interfaceNodeCoord;  //after merging, 

   CRConnectivityPtr   _faceCells;
   CRConnectivityPtr   _faceNodes;
   ArrayVecD3Ptr       _coord;

   ofstream   _debugFile;

   int  _nInterfaceNodes;
   int   _interiorFaceSize;
   MeshList _mesh;   //even it is one-element vector, so compatiable with PartMesh class

};


#endif
