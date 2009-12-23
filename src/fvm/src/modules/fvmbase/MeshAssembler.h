#ifndef _MESHASSEMBLER_H_
#define _MESHASSEMBLER_H_

#include "Mesh.h"
#include "Array.h"
#include "Vector.h"
#include <vector>
#include <set>

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

   void  debug_print();

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

   void setFaceCells();
   void setFaceNodes();

   void countInterfaceNodes();
   
   int  getInnerNodesCount();
   int  getInterfaceNodesDuplicatedCount();
   int  getInterfaceNodesCount();

   
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
   map<int, ArrayIntPtr >   _localNodeToGlobal; //this include ghost cells at interfaces, ghostcells on boundary is set =-1

   set<Mesh::VecD3>  _interfaceNodeCoord;  //after merging, 

   CRConnectivityPtr _faceCells;
   CRConnectivityPtr _faceNodes;

   int  _nInterfaceNodes;

   Mesh* _mesh;

};


#endif
