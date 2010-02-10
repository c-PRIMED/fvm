#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "MeshDismantler.h"
#include "CRConnectivity.h"
#include "MultiField.h"

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif


using namespace std;

MeshDismantler::MeshDismantler( const MeshList& meshList ):_mesh( *meshList.at(0) ), _procID(0)
{
   //assert condition for meshList size
    assert( meshList.size() == 1 ); 
    init();
}

MeshDismantler::~MeshDismantler()
{
   vector< Mesh* >::iterator it_mesh;
   for ( it_mesh = _meshList.begin(); it_mesh != _meshList.end(); it_mesh++)
        delete *it_mesh;
}


void
MeshDismantler::init()
{
#ifdef FVM_PARALLEL
   _procID = MPI::COMM_WORLD.Get_rank();
#endif
   //number of meshes
   _nmesh = _mesh.getNumOfAssembleMesh();
   //giving mesh ids
   int id = 9;
   int dim = _mesh.getDimension();
   //construct meshes
   for ( int n = 0; n < _nmesh; n++ )
      _meshList.push_back( new Mesh( dim, id++) );


   setCellsSite();
   setFacesSite();
   setNodesSite();
   setCellsMapper();
   setFaceCells();
   setNodesMapper();
   setFaceNodes();
   setCoord();
   setMesh();
   debug_print();
}

//setting storagesite for cells
void 
MeshDismantler::setCellsSite()
{
   vector<int> siteGhostCount(_nmesh,0);
   vector<int> siteSelfCount (_nmesh,0);
   //inner cell sweeep
   const StorageSite& cellSite = _mesh.getCells();
   const Array<int>& color     = _mesh.getCellColors();
   for ( int n = 0; n < cellSite.getSelfCount(); n++ )
       siteSelfCount[ color[n] ]++;


   //ghostcells on partition border and boundary cells  sweep
   for ( int n = cellSite.getSelfCount(); n < cellSite.getCount(); n++ )
       siteGhostCount[ color[n] ]++;

   //now find newly emerged ghost cells between meshes  
   //loop over inner faces, if they have two different cell colors, add one ghost cell for each side;
   const CRConnectivity& faceCells = _mesh.getAllFaceCells();
   const StorageSite& faceSite     = _mesh.getInteriorFaceGroup().site;
   for ( int n = 0; n < faceSite.getCount(); n++ ){
         int cell1 = faceCells(n,0);
         int cell2 = faceCells(n,1);
         //check if they are different colors, if so, it is mesh boundary, increment ghostCount for both meshes
         if ( color[ cell1 ] != color[ cell2 ] ){
            siteGhostCount[ color[cell1] ]++;
            siteGhostCount[ color[cell2] ]++;
         }
   }
   //forming cellSites
    for ( int id = 0; id < _nmesh; id++ )
         _cellSite.push_back( StorageSitePtr( new StorageSite(siteSelfCount[id], siteGhostCount[id] ) ) );


}

//setting storagesite for faces
void 
MeshDismantler::setFacesSite()
{
   //loop over all faces, if cells connected to a face has the same color, just add that face to corresponding mesh.
   // if has different colors, that face is counted to add  both sharing meshes
   vector<int> faceCount(_nmesh,0);
   const CRConnectivity& faceCells = _mesh.getAllFaceCells();
   const StorageSite& faceSite     = _mesh.getFaces();
   const Array<int>& color = _mesh.getCellColors();
   for ( int n = 0; n < faceSite.getCount(); n++ ){
       int cell1 = faceCells(n,0);
       int cell2 = faceCells(n,1);
       //check if they are different colors, if so, it is mesh boundary, increment ghostCount for both meshes
       if ( color[ cell1 ] != color[ cell2 ] ){
           faceCount[ color[cell1] ]++;
           faceCount[ color[cell2] ]++;
       } else {
           faceCount[ color[cell2] ]++; // or cell1, cell1 == cell2 in here
       }
   }
   //forming faceSites
   for( int id = 0; id < _nmesh; id++ )
      _faceSite.push_back( StorageSitePtr( new StorageSite(faceCount[id]) ) );


}

//setting Storage site for nodes
void 
MeshDismantler::setNodesSite()
{
   //get nodeCells and look at colors of the cells and incremet nodeCount for each mesh
   //count inner nodes for assembly
   vector<int> nodeCount(_nmesh,0);
   const StorageSite&    nodeSite  = _mesh.getNodes();
   //storing glblNOdeIDs for each mesh
   vector< vector<int> > globalNodeIDs(_nmesh);
   for ( int id = 0; id < _nmesh; id++ )
      globalNodeIDs[id].resize(nodeSite.getCount(),-1);


   const StorageSite&    cellSite  = _mesh.getCells();
   const CRConnectivity& cellNodes = _mesh.getCellNodes();
   const Array<int>& color = _mesh.getCellColors();
   //loop over only inner cells nodes
   for ( int n = 0; n < cellSite.getSelfCount(); n++ ) {
       int nnodes = cellNodes.getCount(n);
       int colorID = color[n];
       for ( int i = 0; i < nnodes; i++ ){
           int glblNodeID = cellNodes(n,i);
           //if it is not visited (=-1)
           if ( globalNodeIDs[colorID][glblNodeID] == -1 ) {
               globalNodeIDs[colorID][glblNodeID] = 1; //(=1) means visited
               nodeCount[colorID]++;
           }
       }
   }
  //pushin in vector field
  for ( int id = 0; id < _nmesh; id++ )
     _nodeSite.push_back( StorageSitePtr( new StorageSite(nodeCount[id]) ) );

}

//gettin localToGlobal and globalToLocal for cell
void
MeshDismantler::setCellsMapper()
{
    //lets create copy cellToGlobal for only inner cells
    for ( int id = 0; id < _nmesh; id++ ){
       const StorageSite& cellSite = *_cellSite.at(id);
       _localCellToGlobal.push_back( ArrayIntPtr( new Array<int>( cellSite.getSelfCount() ) )  );
        Array<int>&  localToGlobal = *_localCellToGlobal[id];
        localToGlobal = -1; //initializer 
    }
    //global to local map ( only inner cells)  
    int cellSelfCount = _mesh.getCells().getSelfCount();
    const Array<int>& color = _mesh.getCellColors();
    _globalCellToLocal.resize ( cellSelfCount, -1);
    _globalCellToMeshID.resize( cellSelfCount, -1); 
    vector<int> localCellCount(_nmesh,0);
    for ( int i = 0;  i < cellSelfCount; i++ ){
        _globalCellToMeshID[i] = color[ i ] ;
        _globalCellToLocal[i] = localCellCount[ color[i] ]++;
    }

}

//getting CRConnectivity faceCells
void
MeshDismantler::setFaceCells()
{
     vector<int> localCellID(_nmesh,0); //track local mesh cell ids
     vector<int> localFaceID(_nmesh,0); //track local mesh face ids

     faceCellsInit( localCellID );
     faceCellsAddInteriorFaces      ( localFaceID );
     faceCellsAddPartitionInterfaces( localFaceID, localCellID );
     faceCellsAddMeshInterfaces     ( localFaceID, localCellID );
     faceCellsAddBoundaryInterfaces ( localFaceID, localCellID );
     faceCellsFinishAdd();

}
//faceCount 
void 
MeshDismantler::faceCellsInit( vector<int>& localCellID )
{    
     //initalize local cellID
     for ( int id = 0; id  < _nmesh; id++ )
         localCellID[id] = _cellSite[id]->getSelfCount();

     //init count and finishCount;
     for ( int id = 0; id < _nmesh; id++ ){
        _faceCells.push_back( CRConnectivityPtr( new CRConnectivity( *_faceSite.at(id), *_cellSite.at(id)) ) );
        _faceCells.at(id)->initCount();

        //addCount, each face share only two cells
        const int cellCount = 2;
        for ( int i = 0; i < _faceSite.at(id)->getCount(); i++ )
           _faceCells.at(id)->addCount(i, cellCount); // face has always two cells
        //finish count
        _faceCells.at(id)->finishCount();
     }
}
//interor face adding
void
MeshDismantler::faceCellsAddInteriorFaces( vector<int>& faceID )
{
     //first add interior faces 
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const FaceGroup& interiorFaceGroup = _mesh.getInteriorFaceGroup();
     const StorageSite& interiorFaceSite = interiorFaceGroup.site;
     for ( int i = 0; i < interiorFaceSite.getCount(); i++ ){
        int cell1 = faceCells(i,0);
        int cell2 = faceCells(i,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 == meshID2 ){ //it means this face interior face
           _faceCells.at( meshID1 )->add( faceID[meshID1], _globalCellToLocal[cell1] );
           _faceCells.at( meshID2 )->add( faceID[meshID2], _globalCellToLocal[cell2] );
           faceID[meshID1]++;
        }
     }
}
//partiion face adding
void
MeshDismantler::faceCellsAddPartitionInterfaces( vector<int>& faceID, vector<int>& localCellID )
{
     //add partition interfaces
     int cellSelfCount  = _mesh.getCells().getSelfCount();
     int interfaceCount = _mesh.getInterfaceGroupCount();
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const FaceGroupList& interfaceGroupList = _mesh.getInterfaceGroups();
     //loop over partition faces (
     for ( int i = 0; i < interfaceCount; i++ ){
         const StorageSite& interiorFaceSite = interfaceGroupList[i]->site;
         int offset = interiorFaceSite.getOffset(); //where to begin face
         int nBeg = offset;
         int nEnd = nBeg + interiorFaceSite.getCount();

         //filling faceOffset and ID
         for ( int id = 0; id < _nmesh; id++ ){
           _interfaceOffset[id].push_back( faceID[id]                );
           _interfaceID    [id].push_back( interfaceGroupList[i]->id );
         }

         for ( int n = nBeg; n < nEnd; n++ ){
             int cell1 = faceCells(n,0);
             int cell2 = faceCells(n,1);
             //if inner cells  take it first, second one ghost cell
             if ( cell1 < cellSelfCount ){
                 int meshID = _globalCellToMeshID[ cell1 ];
                _faceCells.at( meshID )->add( faceID[meshID], _globalCellToLocal[cell1] );
                _faceCells.at( meshID )->add( faceID[meshID], localCellID[meshID]++     );
                 faceID[meshID]++;
             } else {
                 int meshID = _globalCellToMeshID[ cell2 ];
                _faceCells.at( meshID )->add( faceID[meshID], _globalCellToLocal[cell2] );
                _faceCells.at( meshID )->add( faceID[meshID], localCellID[meshID]++     );
                 faceID[meshID]++;
             }
         }

         //filling sizes ( if this interface doesn't involve a mesh, following size will be zero for that mesh
         for ( int id = 0; id < _nmesh; id++ )
            _interfaceSize[id].push_back( faceID[id] - _interfaceOffset[id][i] ); 
 
     }

}
//mesh interface adding
void
MeshDismantler::faceCellsAddMeshInterfaces(vector<int>& faceID, vector<int>& localCellID )
{
     //loop over interfaces to see color difference 
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const FaceGroup& interiorFaceGroup = _mesh.getInteriorFaceGroup();
     const StorageSite& interiorFaceSite = interiorFaceGroup.site;

     vector<int>   countMeshInterface(_nmesh,0);
     for ( int i = 0; i < interiorFaceSite.getCount(); i++ ){
        int cell1 = faceCells(i,0);
        int cell2 = faceCells(i,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 != meshID2 ){ //it means this face  is mesh interface
           countMeshInterface[meshID1]++;
           countMeshInterface[meshID2]++;
        }

     }

      //allocate memory 
     _faceIdentifierList.resize( _nmesh );
     //only sweep interior face to find mesh interface
     for ( int i = 0; i < interiorFaceSite.getCount(); i++ ){
        int cell1 = faceCells(i,0);
        int cell2 = faceCells(i,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 != meshID2 ){ //it means this face  is mesh interface
           _faceIdentifierList[meshID1].insert( pair<int,int>(meshID2,i) );
           _faceIdentifierList[meshID2].insert( pair<int,int>(meshID1,i) );
        }
     }

     //loop over meshes
     for ( int id = 0 ; id < _nmesh ; id++ ){
          const multimap<int,int>& faceIdentifier = _faceIdentifierList[id];
         //loop over all meshinterfaces 
          for ( int key = 0; key < _nmesh; key++ ){ 
               //filling faceOffset, ID and sizes
               int nface = faceIdentifier.count(key);
               if ( nface > 0 ){
                  _interfaceOffset[id].push_back( faceID[id] );
                  _interfaceID    [id].push_back( key        );
                  _interfaceSize  [id].push_back( nface      );
               }

               multimap<int,int>::const_iterator it;
               for ( it = faceIdentifier.equal_range(key).first; it != faceIdentifier.equal_range(key).second; it++ ){
                   int glblFaceID = it->second;
                   int cell1  = faceCells(glblFaceID,0);
                   int cell2  = faceCells(glblFaceID,1);
                   int meshID1 = _globalCellToMeshID[cell1];
                   //int meshID2 = _globalCellToMeshID[cell2];
                   if ( id == meshID1 ){
                     _faceCells.at(id)->add( faceID[id], _globalCellToLocal[cell1]  );
                     _faceCells.at(id)->add( faceID[id],  localCellID[id]++);
                     faceID[id]++;
                   } else {
                     _faceCells.at(id)->add( faceID[id], _globalCellToLocal[cell2]  );
                     _faceCells.at(id)->add( faceID[id],  localCellID[id]++);
                     faceID[id]++;
                   }
               }
          }
     }

}
//boundary interface adding
void 
MeshDismantler::faceCellsAddBoundaryInterfaces( vector<int>& faceID, vector<int>& localCellID )
{
     int cellSelfCount  = _mesh.getCells().getSelfCount();
     int boundaryCount  = _mesh.getBoundaryGroupCount();
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const FaceGroupList& boundaryGroupList = _mesh.getBoundaryFaceGroups();
     //loop over partition faces (
     for ( int i = 0; i < boundaryCount; i++ ){
         const StorageSite& boundaryFaceSite = boundaryGroupList[i]->site;
         int offset = boundaryFaceSite.getOffset(); //where to begin face
         int nBeg = offset;
         int nEnd = nBeg + boundaryFaceSite.getCount();
         //filling faceOffset and ID and Type
         const int bounID = boundaryGroupList[i]->id;
         const string&  bType = boundaryGroupList[i]->groupType;
         for ( int id = 0; id < _nmesh; id++ ){
           _boundaryOffset[id].push_back( faceID[id] );
           _boundaryID    [id].push_back( bounID     );
           _boundaryType  [id].push_back( bType      );
         }
         for ( int n = nBeg; n < nEnd; n++ ){
             int cell1 = faceCells(n,0);
             int cell2 = faceCells(n,1);
             //if inner cells  take it first, second one ghost cell
             if ( cell1 < cellSelfCount ){
                 int meshID = _globalCellToMeshID[ cell1 ];
                _faceCells.at( meshID )->add( faceID[meshID], _globalCellToLocal[cell1] );
                _faceCells.at( meshID )->add( faceID[meshID], localCellID[meshID]++     );
                 faceID[meshID]++;
             } else {
                 int meshID = _globalCellToMeshID[ cell2 ];
                _faceCells.at( meshID )->add( faceID[meshID], _globalCellToLocal[cell2] );
                _faceCells.at( meshID )->add( faceID[meshID], localCellID[meshID]++     );
                 faceID[meshID]++;
            }
         }
         //filling sizes
         for ( int id = 0; id < _nmesh; id++ )
            _boundarySize[id].push_back( faceID[id] - _boundaryOffset[id][i] ); 

     }

}
//finishAdd call
void
MeshDismantler::faceCellsFinishAdd()
{
     //init count and finishCount;
     for ( int id = 0; id < _nmesh; id++ )
        _faceCells.at(id)->finishAdd();
}


//getting CRConnectivity faceNodes
void
MeshDismantler::setFaceNodes()
{
    vector<int> faceID(_nmesh,0); //track local mesh face ids
    faceNodesInit();
    faceNodesAddInteriorFaces      ( faceID );
    faceNodesAddPartitionInterfaces( faceID );
    faceNodesAddMeshInterfaces     ( faceID );
    faceNodesAddBoundaryInterfaces ( faceID );
    faceNodesFinishAdd();
}

//faceNodes finish count
void
MeshDismantler::faceNodesInit()
{
   //init count and finishCount;
   const int nodeCount = _mesh.getAllFaceNodes().getRow()[1] - _mesh.getAllFaceNodes().getRow()[0];
   for ( int id = 0; id < _nmesh; id++ ){
       _faceNodes.push_back( CRConnectivityPtr( new CRConnectivity( *_faceSite.at(id), *_nodeSite.at(id)) ) );
       _faceNodes.at(id)->initCount();

       //addCount, each face share only 
       for ( int i = 0; i < _faceSite.at(id)->getCount(); i++ )
          _faceNodes.at(id)->addCount(i, nodeCount); 
       //finish count
       _faceNodes.at(id)->finishCount();
   }
}
//faceNodes add for interior faces
void
MeshDismantler::faceNodesAddInteriorFaces( vector<int>& faceID )
{
     //first add interior faces 
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const CRConnectivity& faceNodes = _mesh.getAllFaceNodes();
     const FaceGroup& interiorFaceGroup = _mesh.getInteriorFaceGroup();
     const StorageSite& interiorFaceSite = interiorFaceGroup.site;
     const int nodeCount = _mesh.getAllFaceNodes().getRow()[1] - _mesh.getAllFaceNodes().getRow()[0];
     for ( int n = 0; n < interiorFaceSite.getCount(); n++ ){
        int cell1 = faceCells(n,0);
        int cell2 = faceCells(n,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 == meshID2 ){ //it means this face interior face
           for ( int i = 0; i < nodeCount; i++ ){
               int glbNodeID = faceNodes(n,i);
               int lclNodeID = _globalToLocalNodes[glbNodeID][meshID1]; //gives local id 
              _faceNodes.at( meshID1 )->add( faceID[meshID1], lclNodeID );
           }
           faceID[meshID1]++;
        }
     }
}
//faceNodes add for partition faces
void
MeshDismantler::faceNodesAddPartitionInterfaces( vector<int>& faceID )
{
     //add partition interfaces
     int interfaceCount = _mesh.getInterfaceGroupCount();
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const CRConnectivity& faceNodes = _mesh.getAllFaceNodes();
     const FaceGroupList& interfaceGroupList = _mesh.getInterfaceGroups();
     const int nodeCount = _mesh.getAllFaceNodes().getRow()[1] - _mesh.getAllFaceNodes().getRow()[0];
     //loop over partition faces (
     for ( int i = 0; i < interfaceCount; i++ ){
         const StorageSite& interiorFaceSite = interfaceGroupList[i]->site;
         int offset = interiorFaceSite.getOffset(); //where to begin face
         int nBeg = offset;
         int nEnd = nBeg + interiorFaceSite.getCount();
         //loop over faces
         for ( int n = nBeg; n < nEnd; n++ ){
              //loop over nodes 
              int cell1 = faceCells(n,0);
              int meshID = _globalCellToMeshID[ cell1 ];
              for ( int j = 0; j < nodeCount; j++ ){
                  int glbNodeID = faceNodes(n,j);
                  int lclNodeID = _globalToLocalNodes[glbNodeID][meshID]; //gives local id 
                  _faceNodes.at( meshID )->add( faceID[meshID], lclNodeID );
              }
              faceID[meshID]++;
         }
     }

}

void
MeshDismantler::faceNodesAddMeshInterfaces(vector<int>& faceID)
{
     //loop over interfaces to see color difference 
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const CRConnectivity& faceNodes = _mesh.getAllFaceNodes();
     const FaceGroup& interiorFaceGroup = _mesh.getInteriorFaceGroup();
     const StorageSite& interiorFaceSite = interiorFaceGroup.site;
     const int nodeCount = _mesh.getAllFaceNodes().getRow()[1] - _mesh.getAllFaceNodes().getRow()[0];
     //all interior face to search mesh interface
     for ( int n = 0; n < interiorFaceSite.getCount(); n++ ){
        int cell1 = faceCells(n,0);
        int cell2 = faceCells(n,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 != meshID2 ){ //it means this face  is mesh interface
           //face in meshID1
           for ( int i = 0; i < nodeCount; i++ ){
                int glblNodeID = faceNodes(n,i);
                int lclNodeID = _globalToLocalNodes[glblNodeID][meshID1]; //gives local id 
               _faceNodes.at ( meshID1 )->add( faceID[meshID1], lclNodeID );
           }
           faceID[meshID1]++;
           //face in meshID2
           for ( int i = 0; i < nodeCount; i++ ){
                int glblNodeID = faceNodes(n,i);
                int lclNodeID = _globalToLocalNodes[glblNodeID][meshID2]; //gives local id 
               _faceNodes.at ( meshID2 )->add( faceID[meshID2], lclNodeID );
           }
           faceID[meshID2]++;
       }
    }
}


void 
MeshDismantler::faceNodesAddBoundaryInterfaces( vector<int>& faceID )
{
     int boundaryCount  = _mesh.getBoundaryGroupCount();
     const CRConnectivity& faceCells = _mesh.getAllFaceCells(); 
     const CRConnectivity& faceNodes = _mesh.getAllFaceNodes();
     const FaceGroupList& boundaryGroupList = _mesh.getBoundaryFaceGroups();
     const int nodeCount = _mesh.getAllFaceNodes().getRow()[1] - _mesh.getAllFaceNodes().getRow()[0];
     //loop over partition faces (
     for ( int i = 0; i < boundaryCount; i++ ){
         const StorageSite& boundaryFaceSite = boundaryGroupList[i]->site;
         int offset = boundaryFaceSite.getOffset(); //where to begin face
         int nBeg = offset;
         int nEnd = nBeg + boundaryFaceSite.getCount();
         for ( int n = nBeg; n < nEnd; n++ ){
             int cell1 = faceCells(n,0);
             int meshID = _globalCellToMeshID[ cell1 ];
             for ( int j = 0; j < nodeCount; j++ ){
                int glbNodeID = faceNodes(n,j);
                int lclNodeID = _globalToLocalNodes[glbNodeID][meshID]; //gives local id 
                _faceNodes.at( meshID )->add( faceID[meshID], lclNodeID );
             }
             faceID[meshID]++;
         } 
     }
}

void
MeshDismantler::faceNodesFinishAdd()
{
     //init count and finishCount;
     for ( int id = 0; id < _nmesh; id++ )
        _faceNodes.at(id)->finishAdd();
}


//local Nodes to global Nodes 
void
MeshDismantler::setNodesMapper()
{
   //count inner nodes for assembly
   vector<int> nodeCount(_nmesh,0);
   const StorageSite&    nodeSite = _mesh.getNodes();
   //allocate globalToLocaNodes
   _globalToLocalNodes.resize( nodeSite.getCount() );
   //allocating glblNOdeIDs for each mesh
   vector< vector<int> > globalNodeIDs(_nmesh);
   for ( int id = 0; id < _nmesh; id++ )
      globalNodeIDs[id].resize(nodeSite.getCount(),-1);

   const StorageSite&    cellSite  = _mesh.getCells();
   const CRConnectivity& cellNodes = _mesh.getCellNodes();
   const Array<int>& color = _mesh.getCellColors();
   //loop over only inner cells nodes
   for ( int n = 0; n < cellSite.getSelfCount(); n++ ) {
       int nnodes = cellNodes.getCount(n);
       int colorID = color[n];
       for ( int i = 0; i < nnodes; i++ ){
           int glblNodeID = cellNodes(n,i);
           //if it is not visited (=-1)
           if ( globalNodeIDs[colorID][glblNodeID] == -1 ) {
               map<int,int>& nodeMap = _globalToLocalNodes[glblNodeID];
               nodeMap[colorID] = nodeCount[colorID];
               globalNodeIDs[colorID][glblNodeID] = 1;
               nodeCount[colorID]++;
           }
       }
   }

}



//setting coordinates
void
MeshDismantler::setCoord()
{
    //allocation memory for coord
    for ( int id = 0; id < _nmesh; id++ )
       _coord.push_back( ArrayVecD3Ptr(new Array<Mesh::VecD3>(_nodeSite.at(id)->getCount())) );

    const StorageSite& nodeSiteGlbl       = _mesh.getNodes();
    const Array<Mesh::VecD3>&   coordGlbl = _mesh.getNodeCoordinates();
    for ( int i = 0; i < nodeSiteGlbl.getCount(); i++ ){
        const map<int,int>&  colorIDToLocalNode = _globalToLocalNodes[i];
        foreach(const IntMap::value_type& mpos, colorIDToLocalNode){
             int colorID     = mpos.first;
             int localNodeID = mpos.second;
            (*_coord[colorID])[localNodeID] = coordGlbl[i];
        }
    }
}

//setting assembled mesh
void
MeshDismantler::setMesh()
{
    //set sites
    setSites();
    //interior face group
    createInteriorFaceGroup();
    //interface group
    createInterFaceGroup();
    //boundary face group
    createBoundaryFaceGroup();
    //setting coordinates
    createCoords();
    //setting faceNodes CRConnecitivty
    createFaceNodes();
    //setting faceCells CRConnectivity
    createFaceCells();

}

//set StorageSites
void
MeshDismantler::setSites()
{
   for ( int id = 0; id < _nmesh; id++ ){
       StorageSite& faceSite = _meshList.at(id)->getFaces();
       StorageSite& cellSite = _meshList.at(id)->getCells();
       StorageSite& nodeSite = _meshList.at(id)->getNodes();
       //setCounts
       faceSite.setCount( _faceSite.at(id)->getCount() );
       int nGhost = _cellSite.at(id)->getCount()-_cellSite.at(id)->getSelfCount();
       cellSite.setCount( _cellSite.at(id)->getSelfCount(), nGhost );
       nodeSite.setCount( _nodeSite.at(id)->getCount() );
   }
}
//interior face
void
MeshDismantler::createInteriorFaceGroup()
{
    for ( int id = 0; id < _nmesh; id++ ){
         int nInteriorFace = 0;
         const CRConnectivity& faceCells = *_faceCells.at(id);
         const StorageSite&    faceSite  = *_faceSite.at(id);
         const int cellSelfCount  = _cellSite.at(id)->getSelfCount();
         for ( int i = 0; i < faceSite.getCount(); i++ ){
              int cell1  = faceCells(i,0); 
              int cell2  = faceCells(i,1);
              if ( (cell1 < cellSelfCount) && (cell2 < cellSelfCount) )
                 nInteriorFace++;
         }
        _meshList.at(id)->createInteriorFaceGroup( nInteriorFace );
    }

}
//interface group
void
MeshDismantler::createInterFaceGroup()
{
    for ( int id = 0; id < _nmesh; id++ ){
         for ( unsigned int i = 0; i < _interfaceOffset[id].size(); i++ ){
             const int size   = _interfaceSize[id][i];
             const int offset = _interfaceOffset[id][i];
             const int interfaceID     = _interfaceID[id][i];
             if ( size > 0 )
                _meshList.at(id)->createInterfaceGroup( size, offset, interfaceID);
         }
     }
}
//boundary face group
void  
MeshDismantler::createBoundaryFaceGroup()
{
    for ( int id = 0; id < _nmesh; id++ ){
        for ( unsigned int i = 0; i < _boundaryOffset[id].size(); i++ ){
            const int size   = _boundarySize[id][i];
            const int offset = _boundaryOffset[id][i];
            const int boundaryID     = _boundaryID[id][i];
            const string& bType = _boundaryType[id][i];
            if ( size > 0 )
                _meshList.at(id)->createBoundaryFaceGroup( size, offset, boundaryID, bType );
        }
    }
}

//setting coords
void
MeshDismantler::createCoords()
{
   for ( int id = 0 ; id < _nmesh; id++ )
      _meshList.at(id)->setCoordinates( _coord.at(id) );
}
//set faceNodes for meshList
void
MeshDismantler::createFaceNodes()
{
    for ( int id = 0; id < _nmesh; id++ )
        _meshList.at(id)->setFaceNodes( _faceNodes.at(id) );
}
//set faceCells for meshList
void
MeshDismantler::createFaceCells()
{
    for ( int id = 0; id < _nmesh; id++ )
        _meshList.at(id)->setFaceCells( _faceCells.at(id) );
}


void  
MeshDismantler::debug_print()
{
   debug_cell_site();
   debug_face_site();
   debug_node_site();
   debug_cells_mapper();
   debug_face_cells();
   debug_nodes_mapper();
   debug_face_nodes();
}


void
MeshDismantler::debug_cell_site()
{
    debug_file_open("cellSite");
    for ( int id = 0; id < _nmesh; id++ )
        _debugFile <<"meshid = " << id <<  "   selfCount = " << _cellSite.at(id)->getSelfCount() << "   count = " << _cellSite.at(id)->getCount() << endl;
    debug_file_close();
}

void
MeshDismantler::debug_face_site()
{
    debug_file_open("faceSite");
    for ( int id = 0; id < _nmesh; id++ )
        _debugFile <<"meshid = " << id <<  "   count = " << _faceSite.at(id)->getCount() << endl;
    debug_file_close();
}

void
MeshDismantler::debug_node_site()
{
    debug_file_open("nodeSite");
    for ( int id = 0; id < _nmesh; id++ )
        _debugFile <<"meshid = " << id <<  "   count = " << _nodeSite.at(id)->getCount() << endl;
    debug_file_close();
}

void
MeshDismantler::debug_cells_mapper()
{ 
    debug_file_open("cellsMapper");
    for ( unsigned int i = 0; i < _globalCellToMeshID.size(); i++ )
        _debugFile << "glblID = " << i << "   meshID  = " << _globalCellToMeshID[i] << endl;
    _debugFile << endl;
    for ( unsigned int i = 0; i < _globalCellToLocal.size(); i++ )
        _debugFile << "glblID = " << i << "   localID = " << _globalCellToLocal[i] << endl;
    debug_file_close();

}

void
MeshDismantler::debug_nodes_mapper()
{ 
    debug_file_open("nodesMapper");
    for ( unsigned i = 0; i < _globalToLocalNodes.size(); i++ ){
       const map<int,int>& nodeMap = _globalToLocalNodes[i];
       foreach(const IntMap::value_type& mpos, nodeMap){
           int colorID     = mpos.first;
           int localNodeID = mpos.second;
          _debugFile << "glblNodeID = " << i << "   meshID = " <<  colorID << "   localNodeID = " << localNodeID << endl;
       }
    }
    debug_file_close();
}

void
MeshDismantler::debug_face_cells()
{
    debug_file_open("faceCells");
    for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite   & faceSite  = *_faceSite.at(id);
        const CRConnectivity& faceCells = *_faceCells.at(id);
        _debugFile << " meshID : " << id << endl;
        for ( int n = 0; n < faceSite.getCount(); n++ ){
            _debugFile << "faceCells(" << n << " ) = ";
            for ( int i = 0; i < faceCells.getCount(n); i++ ){
                _debugFile << faceCells(n,i) << "     ";
            }
            _debugFile << endl;
        }
    }
    debug_file_close();

}

void
MeshDismantler::debug_face_nodes()
{
    debug_file_open("faceNodes");
    for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite   & faceSite  = *_faceSite.at(id);
        const CRConnectivity& faceNodes = *_faceNodes.at(id);
        _debugFile << " meshID : " << id << endl;
        for ( int n = 0; n < faceSite.getCount(); n++ ){
            _debugFile << "faceNodes(" << n << " ) = ";
            for ( int i = 0; i < faceNodes.getCount(n); i++ ){
                _debugFile << faceNodes(n,i) << "     ";
            }
            _debugFile << endl;
        }
    }
    debug_file_close();

}


void
MeshDismantler::debug_file_open( const string& fname_ )
{  
     stringstream ss;
     ss << _procID;

     string  fname = "MESHDISMANTLER_"+fname_+"_proc"+ss.str()+".dat";
    _debugFile.open( fname.c_str() );
}

void
MeshDismantler::debug_file_close()
{
    _debugFile.close();
}