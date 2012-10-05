// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include "MeshDismantler.h"
#include "CRConnectivity.h"
#include "MultiField.h"




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
   
  _cellFaces = _mesh.getAllFaceCells().getTranspose();
#ifdef FVM_PARALLEL
   _procID = MPI::COMM_WORLD.Get_rank();
   _nPart = MPI::COMM_WORLD.Get_size();
#endif
   //number of meshes
   _nmesh = _mesh.getNumOfAssembleMesh();
   //giving mesh ids
   int dim = _mesh.getDimension();
   //construct meshes	
   for ( int n = 0; n < _nmesh; n++ )
      _meshList.push_back( new Mesh( dim) );


   setCellsSite();
   setFacesSite();
   setNodesSite();
   setSites();
   setCellsMapper();
   setFaceCells();
   setNodesMapper();
   setFaceNodes();
   setCoord();
   setMesh();
   setMappers();
   set_local_global();

   for ( int n = 0; n < _nmesh; n++ ){
      const StorageSite& cellSite = _meshList.at(n)->getCells();
      const StorageSite& faceSite = _meshList.at(n)->getFaces();
      _meshList.at(n)->eraseConnectivity(cellSite, cellSite);
      _meshList.at(n)->eraseConnectivity(cellSite, faceSite);
      //uniquie
      _meshList.at(n)->uniqueFaceCells();
   }  

   
   setCellCellsGhostExt();


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

    //we will first identify newly emerged cellcell2 cells
    vector< set<int> >  gatherCells(_nmesh);
    for ( int n = 0; n < faceSite.getCount(); n++ ){
       int cell1 = faceCells(n,0);
       int cell2 = faceCells(n,1);
       //check if they are different colors, if so, it is mesh boundary, increment ghostCount for both meshes
       if ( color[ cell1 ] != color[ cell2 ] ){
          gatherCells[ color[cell1] ].insert(cell2);
          gatherCells[ color[cell2] ].insert(cell1);
       }
    }
    //now loop over gather cells and check its global cellCells connectivity
    const multimap<int,int>& cellCellsGlobal = _mesh.getCellCellsGlobal();
    const map<int,int>&      globalToLocal   = _mesh.getGlobalToLocal();
    for ( int id = 0; id < _nmesh; id++ ){
      const set<int>& cells = gatherCells[id];
      set<int> cells2; //storing cellCells2 cells
      int countLevel1 = _cellSite[id]->getCount();
      //loop over gather cells
      foreach(const set<int>::value_type& mpos, cells){
          int cellID = mpos;
	  multimap<int,int>::const_iterator it;
          for ( it = cellCellsGlobal.equal_range(cellID).first; it != cellCellsGlobal.equal_range(cellID).second; it++ ){
	     const int localID = globalToLocal.find(it->second)->second;
	     //if it is not gathercells, we accep it, and color[localID] make sure that it doesn't pick inner cells 
	     //and cells2 make sure that this is not counted twice
	     if ( cells.count(localID) == 0 &&  color[localID] != id && cells2.count(localID) == 0){
	        countLevel1++;
		cells2.insert(localID);
             }
          }
      }
      //update countLevel1
      // _cellSite[id]->setCountLevel1(countLevel1);
   }
    

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
     faceCellsAddBoundaryInterfaces ( localFaceID, localCellID );
     faceCellsAddMeshInterfaces     ( localFaceID, localCellID );
     faceCellsAddPartitionInterfaces( localFaceID, localCellID );
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
         const StorageSite& faceSite = _meshList.at(id)->getFaces();
         const StorageSite& cellSite = _meshList.at(id)->getCells();
        _faceCells.push_back( CRConnectivityPtr( new CRConnectivity( faceSite, cellSite ) ) );
        _faceCells.at(id)->initCount();

        //addCount, each face share only two cells
        const int cellCount = 2;
        for ( int i = 0; i < faceSite.getCount(); i++ )
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
           _interfaceID    [id].push_back( -interfaceGroupList[i]->id );
         }

         for ( int n = nBeg; n < nEnd; n++ ){
             int cell1 = faceCells(n,0);
             int cell2 = faceCells(n,1);
             //if inner cells  take it first, second one ghost cell
             if ( cell1 < cellSelfCount ){
                 int meshID = _globalCellToMeshID[ cell1 ];
                _faceCells.at( meshID )->add( faceID[meshID], _globalCellToLocal[cell1] );
                _faceCells.at( meshID )->add( faceID[meshID], localCellID[meshID]       );
                 _globalToLocalFaces[meshID][n] = faceID[meshID];
                 localCellID[meshID]++;
                 faceID[meshID]++;
             } else {
                 int meshID = _globalCellToMeshID[ cell2 ];
                _faceCells.at( meshID )->add( faceID[meshID], _globalCellToLocal[cell2] );
                _faceCells.at( meshID )->add( faceID[meshID], localCellID[meshID]     );
                 _globalToLocalFaces[meshID][n] = faceID[meshID];
                 localCellID[meshID]++;
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
//                   int meshID2 = _globalCellToMeshID[cell2];
                   if ( id == meshID1 ){
                     _faceCells.at(id)->add( faceID[id], _globalCellToLocal[cell1]  );
                     _faceCells.at(id)->add( faceID[id],  localCellID[id]);
                     _globalToLocalFaces[id][glblFaceID] = faceID[id];
                     localCellID[id]++;
                     faceID[id]++;
                   } else {
                     _faceCells.at(id)->add( faceID[id], _globalCellToLocal[cell2]  );
                     _faceCells.at(id)->add( faceID[id],  localCellID[id] );
                     _globalToLocalFaces[id][glblFaceID] = faceID[id];
                     localCellID[id]++;
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
    _globalToLocalFaces.resize( _nmesh );
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
    faceNodesAddBoundaryInterfaces ( faceID );
    faceNodesAddMeshInterfaces     ( faceID );
    faceNodesAddPartitionInterfaces( faceID );
    faceNodesFinishAdd();
}

//faceNodes finish count
void
MeshDismantler::faceNodesInit()
{
   //init count and finishCount;
   for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite&  faceSite = _meshList.at(id)->getFaces();
        const StorageSite&  nodeSite = _meshList.at(id)->getNodes();
       _faceNodes.push_back( CRConnectivityPtr( new CRConnectivity( faceSite, nodeSite ) ) );
       _faceNodes.at(id)->initCount();
   }


     //first add interior faces addCounts
     vector<int> faceIndx(_nmesh,0);
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const CRConnectivity& faceNodes = _mesh.getAllFaceNodes();
     const FaceGroup& interiorFaceGroup = _mesh.getInteriorFaceGroup();
     const StorageSite& interiorFaceSite = interiorFaceGroup.site;
     for ( int n = 0; n < interiorFaceSite.getCount(); n++ ){
        int cell1 = faceCells(n,0);
        int cell2 = faceCells(n,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 == meshID2 ){ //it means this face interior face
            const int nodeCount = faceNodes.getCount(n);
            _faceNodes.at( meshID1 )->addCount( faceIndx[meshID1]++, nodeCount );
        }
     }

     //add partition interfaces to addCount
     const int interfaceCount = _mesh.getInterfaceGroupCount();
     const FaceGroupList& interfaceGroupList = _mesh.getInterfaceGroups();
     //loop over partition faces (
     for ( int i = 0; i < interfaceCount; i++ ){
         const StorageSite& interiorFaceSite = interfaceGroupList[i]->site;
         int offset = interiorFaceSite.getOffset(); //where to begin face
         int nBeg = offset;
         int nEnd = nBeg + interiorFaceSite.getCount();
         //loop over faces
         for ( int n = nBeg; n < nEnd; n++ ){
              //loop over nodes 
              const int cell1 = faceCells(n,0);
              const int meshID = _globalCellToMeshID[ cell1 ];
              const int nodeCount = faceNodes.getCount(n);
             _faceNodes.at( meshID )->addCount( faceIndx[meshID]++, nodeCount );
         }
     }

      //loop over mesh interfaces to addCounts
     //all interior face to search mesh interface
     for ( int n = 0; n < interiorFaceSite.getCount(); n++ ){
        int cell1 = faceCells(n,0);
        int cell2 = faceCells(n,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 != meshID2 ){ //it means this face  is mesh interface
           //face in meshID1
           const int nodeCount = faceCells.getCount(n);
           _faceNodes.at ( meshID1 )->addCount( faceIndx[meshID1]++, nodeCount ); 
           //face in meshID2
           _faceNodes.at ( meshID2 )->addCount( faceIndx[meshID2]++, nodeCount ); 
         }
      }

     //loop over boundary faces to addCount
     const int boundaryCount  = _mesh.getBoundaryGroupCount();
     const FaceGroupList& boundaryGroupList = _mesh.getBoundaryFaceGroups();
     //loop  boundary faces to addCount
     for ( int i = 0; i < boundaryCount; i++ ){
         const StorageSite& boundaryFaceSite = boundaryGroupList[i]->site;
         const int offset = boundaryFaceSite.getOffset(); //where to begin face
         const int nBeg = offset;
         const int nEnd = nBeg + boundaryFaceSite.getCount();
         for ( int n = nBeg; n < nEnd; n++ ){
             const int cell1 = faceCells(n,0);
             const int meshID = _globalCellToMeshID[ cell1 ];
             const int nodeCount = faceNodes.getCount(n);
             _faceNodes.at( meshID )->addCount( faceIndx[meshID]++, nodeCount );
         } 
     }


    //finish count
    for ( int id = 0; id < _nmesh; id++ )
       _faceNodes.at(id)->finishCount();
 
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
     for ( int n = 0; n < interiorFaceSite.getCount(); n++ ){
        int cell1 = faceCells(n,0);
        int cell2 = faceCells(n,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 == meshID2 ){ //it means this face interior face
           const int nodeCount = faceCells.getCount(n);
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
              const int nodeCount = faceCells.getCount(n);
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
     //all interior face to search mesh interface
     for ( int n = 0; n < interiorFaceSite.getCount(); n++ ){
        int cell1 = faceCells(n,0);
        int cell2 = faceCells(n,1);
        int meshID1 = _globalCellToMeshID[cell1];
        int meshID2 = _globalCellToMeshID[cell2];
        if ( meshID1 != meshID2 ){ //it means this face  is mesh interface
           const int nodeCount = faceCells.getCount(n);
           //face in meshID1
           for ( int i = 0; i < nodeCount; i++ ){
                int glblNodeID = faceNodes(n,i);
                int lclNodeID = _globalToLocalNodes[glblNodeID][meshID1]; //gives local id 
               _faceNodes.at ( meshID1 )->add( faceID[meshID1], lclNodeID );
           }
           //face in meshID2 (consistent with global faceNodes connectivity)
            for ( int i = nodeCount-1; i >= 0; i-- ){
                int glblNodeID = faceNodes(n,i);
                int lclNodeID = _globalToLocalNodes[glblNodeID][meshID2]; //gives local id 
               _faceNodes.at ( meshID2 )->add( faceID[meshID2], lclNodeID );
           }
           faceID[meshID1]++;
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
     //loop over partition faces (
     for ( int i = 0; i < boundaryCount; i++ ){
         const StorageSite& boundaryFaceSite = boundaryGroupList[i]->site;
         int offset = boundaryFaceSite.getOffset(); //where to begin face
         int nBeg = offset;
         int nEnd = nBeg + boundaryFaceSite.getCount();
         for ( int n = nBeg; n < nEnd; n++ ){
             int cell1 = faceCells(n,0);
             int meshID = _globalCellToMeshID[ cell1 ];
              const int nodeCount = faceCells.getCount(n);
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
    //interior face group
    createInteriorFaceGroup();
    //boundary face group
    createBoundaryFaceGroup();
    //interface group
    createInterFaceGroup();
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
       cellSite.setCountLevel1( _cellSite.at(id)->getCountLevel1() );
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
//fill scatter and gather maps
void
MeshDismantler::setMappers()
{
       partitionInterfaceMappers();
       meshInterfaceMappers();

}

//partition interface mappers
void
MeshDismantler::partitionInterfaceMappers()
{
     //mappers
     const StorageSite::ScatterMap& scatterMap = _mesh.getCells().getScatterMap();
     const StorageSite::GatherMap&  gatherMap  = _mesh.getCells().getGatherMap ();
     //interfaceList
     const FaceGroupList& faceGroupList = _mesh.getInterfaceGroups();
     //interfacecount
     int interfaceGroupCount = _mesh.getInterfaceGroupCount();
     //loop over interfaces (only partition interfaces)
     for ( int i = 0; i < interfaceGroupCount; i++ ) {
         //get key for f
         int interfaceID = -faceGroupList[i]->id;
         //mappers in Storagesite is held by the same storage site, so this site 
         //can be used to in scatterMap and gatherMap
         const StorageSite& interfaceSite = faceGroupList[i]->site;
         //from intrface ID and meshID (since this is from meshID), we can get that specific storagesite
         const  int meshID0 = 0;
         Mesh::PartIDMeshIDPair  pairID = make_pair<int,int>( interfaceID, meshID0 );
         const StorageSite* ghostSite = _mesh.getGhostCellSiteScatter( pairID );
         //get mapper arrays
         const Array<int>& scatterArray = *scatterMap.find(ghostSite)->second;
         const Array<int>& gatherArray  = *gatherMap.find(ghostSite)->second;
         EntryVecMap   scatterArrayMap;
         EntryVecMap   gatherArrayMap;
         getGatherArrays ( gatherArray , gatherArrayMap , interfaceSite );
         getScatterArrays( scatterArray, scatterArrayMap, interfaceSite );

         //get scatterArrayMap
         foreach ( const EntryVecMap::value_type& pos, scatterArrayMap ){
              const pair<int,int>& entry = pos.first;
              const int meshID = entry.first;
              const int otherMeshID  = entry.second;
              const int size   = int(pos.second.size());
              //reference ot mappers to fill in
              StorageSite::ScatterMap& scatterMapLocal = _meshList.at(meshID)->getCells().getScatterMap();
              StorageSite::GatherMap& gatherMapLocal   = _meshList.at(meshID)->getCells().getGatherMap();
              //storagesite (used for both scatter and gathersites), pass parent (getCells())
              shared_ptr<StorageSite> siteScatterLocal( new StorageSite(size) );

              //copy scatterArray to Array<int> 
              int scatterSize =  int(scatterArrayMap[entry].size());
              ArrayIntPtr scatterArrayLocal( new Array<int>( scatterSize ) );
              for ( int i = 0; i < scatterSize; i++ ) 
                   (*scatterArrayLocal)[i] = scatterArrayMap[entry][i];
              //copy gatherArray to Array<int>
              int gatherSize =  int(gatherArrayMap[entry].size());
              ArrayIntPtr gatherArrayLocal( new Array<int>( gatherSize ) );
              for ( int i = 0; i < gatherSize; i++ ) 
                   (*gatherArrayLocal)[i] = gatherArrayMap[entry][i];

              //setting scatterID and gatherID
              siteScatterLocal->setScatterProcID( ghostSite->getScatterProcID() );
              siteScatterLocal->setGatherProcID ( ghostSite->getGatherProcID()  );
              //setting tag (shifting 16 bits to left)
              int packed_info = (std::max(meshID,otherMeshID) << 16 ) | ( std::min(meshID,otherMeshID) );
              siteScatterLocal->setTag( packed_info );	

              Mesh::PartIDMeshIDPair  pairID = make_pair<int,int>( ghostSite->getGatherProcID(), otherMeshID);
              _meshList.at(meshID)->createGhostCellSiteScatter( pairID, siteScatterLocal );
              _meshList.at(meshID)->createGhostCellSiteGather ( pairID, siteScatterLocal );
              scatterMapLocal[ siteScatterLocal.get() ] = scatterArrayLocal;
              gatherMapLocal [ siteScatterLocal.get() ] = gatherArrayLocal;
         }
      }


}


void
MeshDismantler::getScatterArrays( const Array<int>& scatterArray, EntryVecMap& scatterArrayLocal,  const StorageSite& site )
{
    // get counts for each mesh
     EntryMap sizeScatter;
     int offset = site.getOffset(); //where to begin face
     int iBeg = offset;
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const Array<int>&  color        = _mesh.getCellColors();
     const Array<int>&  colorOther   = _mesh.getCellColorsOther();
     //loop over partition faces
     for ( int i = 0; i < scatterArray.getLength(); i++ ){
         const int faceIndx = iBeg+i;
         //check gather cellColor
         const int cell1  = faceCells(faceIndx,1); //outside cell
         const int thisMeshID  = color[ scatterArray[i] ];
         const int otherMeshID = colorOther[ cell1 ];
         EntryIndex eIndex = make_pair<int,int>(thisMeshID, otherMeshID);
         //sizeScatter[eIndex]++;
         const int  cellID = _globalCellToLocal [ scatterArray[i] ];

         scatterArrayLocal[eIndex].push_back(cellID);
     }

}


void
MeshDismantler::getGatherArrays(const Array<int>& gatherArray, EntryVecMap& gatherArrayLocal, const StorageSite& site )
{
      // get counts for each mesh
     EntryMap sizeScatter;
     int offset = site.getOffset(); //where to begin face
     int iBeg = offset;
     const CRConnectivity& faceCells = _mesh.getAllFaceCells();
     const Array<int>&  color        = _mesh.getCellColors();
     const Array<int>&  colorOther   = _mesh.getCellColorsOther();
     //loop over partition faces
     for ( int i = 0; i < gatherArray.getLength(); i++ ){
         const int faceIndx = iBeg+i;
         //check gather cellColor
         const int cell0  = faceCells(faceIndx,0); //outside cell
         const int thisMeshID  = color[  cell0 ];
         const int otherMeshID = colorOther[ gatherArray[i] ];
         EntryIndex eIndex = make_pair<int,int>(thisMeshID, otherMeshID);
         const int glblFaceID  = (*_cellFaces)(gatherArray[i],0);
         const int localFaceID = _globalToLocalFaces[thisMeshID][glblFaceID];
         const int  cellID = (*_faceCells[thisMeshID])(localFaceID,1); // 0 inner , 1 ghost cells        
         gatherArrayLocal[eIndex].push_back(cellID);
     }

}


//mesh interfaces
void
MeshDismantler::meshInterfaceMappers()
{
     //loop over meshes
     for ( int id = 0 ; id < _nmesh ; id++ ){
          const multimap<int,int>& faceIdentifier  = _faceIdentifierList[id];
          StorageSite::GatherMap & gatherMapLocal  = _meshList.at(id)->getCells().getGatherMap();
          const StorageSite* thisSite = &_meshList.at(id)->getCells();
          //loop over all meshinterfaces  (key = all other meshes)
          for ( int key = 0; key < _nmesh; key++ ){ 
               //filling scatter (on this mesh) and gather Array (on other mesh)
               int nface = faceIdentifier.count(key);
               if ( nface > 0 ){
                   StorageSite::ScatterMap& scatterMapLocal = _meshList.at(key)->getCells().getScatterMap(); 
                   const StorageSite* otherSite = &_meshList.at(key)->getCells();
                   ArrayIntPtr scatterArrayLocal( new Array<int>( nface ) ); //for other side mesh
                   ArrayIntPtr gatherArrayLocal ( new Array<int>( nface ) ); //for this side mesh
                   multimap<int,int>::const_iterator it;
                   int indx = 0;
                   for ( it = faceIdentifier.equal_range(key).first; it != faceIdentifier.equal_range(key).second; it++ ){
                       //fill this mesh gather
                       const int glblFaceID  = it->second;
                       int localFaceID = _globalToLocalFaces[id][glblFaceID];
                       const int gatherCellID = (*_faceCells[id])(localFaceID,1);
                       (*gatherArrayLocal)[indx] = gatherCellID;
                       //now other mesh to fill scatter arrays
                       localFaceID = _globalToLocalFaces[key][glblFaceID];
                       const int scatterCellID = (*_faceCells[key])(localFaceID,0);
                       (*scatterArrayLocal)[indx] = scatterCellID;
                       indx++;
                   }
                   gatherMapLocal.insert ( make_pair(otherSite, gatherArrayLocal ) ); //gather this side mesh, so we key with other site StorageSite*
                   scatterMapLocal.insert( make_pair(thisSite , scatterArrayLocal) ); //scatter other side mesh, so we key with this site StorageSite*
               }
          }
     }


}

//filling mesh localToGlobal and globalToLocal
void 
MeshDismantler::set_local_global()
{
    const int nmesh = int( _meshList.size() );
    //creating cellID MultiField to use sync() operation
    shared_ptr<MultiField>  cellMultiField = shared_ptr<MultiField>( new MultiField()    );
    shared_ptr<Field>       cellField      = shared_ptr<Field>     ( new Field("globalcellID") );
 
    for ( int id = 0; id < nmesh; id++ ){
       const StorageSite* site = &_meshList[id]->getCells();
       MultiField::ArrayIndex ai( cellField.get(), site );
       shared_ptr<Array<int> > cIndex(new Array<int>(site->getCountLevel1()));
       *cIndex = -1;
       cellMultiField->addArray(ai,cIndex);
    }

    //global numbering 
    const int globalOffset = global_offset();
    int offset = globalOffset;
    for ( int id = 0; id < nmesh; id++ ){
       const Mesh&    mesh = *_meshList.at(id);
       const StorageSite* site = &_meshList[id]->getCells();
       MultiField::ArrayIndex ai( cellField.get(), site );
       Array<int>&  localCell = dynamic_cast< Array<int>& >( (*cellMultiField)[ai] ); 
       //global numbering inner cells
       const int selfCount = site->getSelfCount();
       for ( int i = 0; i < selfCount; i++ )
          localCell[i] = offset + i;
       //update offset 
       offset += selfCount;
       //loop over boundaries and global number boundary cells
       const FaceGroupList&  bounGroupList = mesh.getBoundaryFaceGroups();
       const CRConnectivity& faceCells  = mesh.getAllFaceCells();
       for ( int i = 0; i < mesh.getBoundaryGroupCount(); i++ ){
          const int ibeg = bounGroupList[i]->site.getOffset();
          const int iend = ibeg + bounGroupList[i]->site.getCount();
          int indx=0;
          for ( int ii = ibeg; ii < iend; ii++ ){
             localCell[ faceCells(ii,1)] = offset + indx;
             indx++;
          }
          //update offset
          offset += iend-ibeg;
       }
       
    }

    //sync opeartion
    cellMultiField->sync();

    //create localToGlobal array and assign it in Mesh
    for ( int id = 0; id < nmesh; id++ ){
       Mesh& mesh = *_meshList.at(id);
       mesh.createLocalGlobalArray();
       const StorageSite* site = &_meshList[id]->getCells();
       MultiField::ArrayIndex ai( cellField.get(), site );
       const Array<int>&  localCell = dynamic_cast< const Array<int>& >( (*cellMultiField)[ai] ); 
       Array<int>& localToGlobal = mesh.getLocalToGlobal();
       for ( int i = 0; i < localCell.getLength(); i++ ){
          localToGlobal[i] = localCell[i];
          assert( localCell[i] != -1 );
       }

       //copying GlobalToLocal
       map<int,int>& globalToLocal = mesh.getGlobalToLocal();
       for ( int i = 0; i < localCell.getLength(); i++ ){
          globalToLocal[ localToGlobal[i] ] = i;
       }
    }

}



//get offset value for global numbering for each partition
int
MeshDismantler::global_offset()
{
   const int nmesh = int( _meshList.size() );
   int count = 0;
   //get offsets for 
   for ( int id = 0; id < nmesh; id++ ){
      const Mesh&    mesh = *_meshList.at(id);
      const StorageSite& cellSite = mesh.getCells();
      const FaceGroupList&  bounGroupList = mesh.getBoundaryFaceGroups();
      int bounCount = 0;
      for ( int i = 0; i < mesh.getBoundaryGroupCount(); i++ )
         bounCount += bounGroupList[i]->site.getCount();
      const int selfCount = cellSite.getSelfCount();
      count += selfCount + bounCount;
   }
   
   
    //allocation holding each partiton offset
    int *counts = new int[ _nPart ];
   //MPI calls allgather to know offsets
#ifdef FVM_PARALLEL   
   MPI::COMM_WORLD.Allgather( &count, 1, MPI::INT, counts, 1, MPI::INT);
#endif   
   
   //compute offsets for each partition
   int offset = 0;
   for ( int i = 0; i < _procID; i++ )
      offset += counts[i];

   //delete allocation counts
   delete [] counts;
   return offset;

}


   
void
MeshDismantler::setCellCellsGhostExt()
{
    //create sendCounts
    for (int n=0; n<_nmesh; n++ ){
       Mesh& mesh = *_meshList[n];
       mesh.createScatterGatherCountsBuffer();
       //partitioner interfaces
       mesh.syncCounts();
     }  

    for (int n=0; n<_nmesh; n++ ){
       Mesh& mesh = *_meshList[n];
       //mesh interfaces
       mesh.recvScatterGatherCountsBufferLocal();
    }
    
    
    //create indices
    for (int n=0; n<_nmesh; n++ ){
       Mesh& mesh = *_meshList[n];
       mesh.createScatterGatherIndicesBuffer();
       //partitioner interfaces
       mesh.syncIndices();
    }       
    
    for (int n=0; n<_nmesh; n++ ){
       Mesh& mesh = *_meshList[n];
       //mesh interfaces
       mesh.recvScatterGatherIndicesBufferLocal();
    }

    //creaet cellCellsGhostExt
    for (int n=0; n<_nmesh; n++ ){
       Mesh& mesh = *_meshList[n];
       //mesh interfaces
       mesh.createCellCellsGhostExt();
    }
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
   debug_scatter_mappers();
   debug_gather_mappers();
}


void
MeshDismantler::debug_cell_site()
{
    debug_file_open("cellSite");
    for ( int id = 0; id < _nmesh; id++ )
        _debugFile <<"meshid = " << id <<  "   selfCount = " << _cellSite.at(id)->getSelfCount() << "   count = " << _cellSite.at(id)->getCount() <<
                     "   countLevel1 = " << _cellSite.at(id)->getCountLevel1() <<endl;
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
MeshDismantler::debug_scatter_mappers()
{
    debug_file_open("scatterMappers");
    //creating mappers between cell storage site to mesh id
    map< const StorageSite*, int > siteMeshMapper; //key  = storage site, value = mesh ID of cellSite
    for ( int id = 0; id < _nmesh; id++ )
        siteMeshMapper[&_meshList.at(id)->getCells()] = id;

    //first make sure consistent order for scatterMap since map order with key which is pointer
    //first mesh interface (increasing order) then partition interface (increasing order)
    vector< const StorageSite* > scatterSiteVec;
    map <int, int> packIDToMeshID;
    vector< int > scatterSiteVecID;
    IntStorageSiteMap   orderedMeshInterface;
    for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite::ScatterMap&  scatterMap = _meshList.at(id)->getCells().getScatterMap();
        foreach ( const StorageSite::ScatterMap::value_type& pos, scatterMap ){
            const StorageSite& scatterSite = *pos.first;
            if ( scatterSite.getGatherProcID() == -1 ){  //this means mesh interface
                int neighMeshID = siteMeshMapper[ pos.first ];
                const int packed_info = ( (id  << 16)  | neighMeshID );
                orderedMeshInterface.insert( make_pair<int, const StorageSite*>(packed_info, pos.first) );
                packIDToMeshID.insert( make_pair<int,int>(packed_info,id) );
            }
        }
    }

   //fill scatterSiteVec, map data structure already ordered them
   foreach ( const IntStorageSiteMap::value_type& pos, orderedMeshInterface ){
          scatterSiteVec.push_back( pos.second );
          const int meshID = packIDToMeshID[pos.first];
          scatterSiteVecID.push_back( meshID );
  }
   //clear for usage in partition interface
    packIDToMeshID.clear();
    //order partition interface (map will reorder)
    IntStorageSiteMap orderedPartInterface;
    for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite::ScatterMap&  scatterMap = _meshList.at(id)->getCells().getScatterMap();
        foreach ( const StorageSite::ScatterMap::value_type& pos, scatterMap ){
            const StorageSite& scatterSite = *pos.first;
            if ( scatterSite.getGatherProcID() != -1 ){  //this means partition face
                const int neighID = scatterSite.getGatherProcID();
                const int tag     = scatterSite.getTag();
                const int packed_info = ( (id << 31)  | (neighID  << 24) | tag );
                orderedPartInterface.insert( make_pair<int, const StorageSite*>(packed_info, pos.first) );
                packIDToMeshID.insert( make_pair<int,int>(packed_info,id) );
            }
        }
     }

    //fill scatterSiteVec, map data structure already 
    foreach ( const IntStorageSiteMap::value_type& pos, orderedPartInterface ){
          scatterSiteVec.push_back( pos.second );
          const int meshID = packIDToMeshID[pos.first];
          scatterSiteVecID.push_back( meshID );
    }

     int indx = 0;
     foreach( const vector<const StorageSite*>::value_type& pos, scatterSiteVec){
          if ( pos->getGatherProcID() != -1 ){  //this means partition face
               const int id = scatterSiteVecID[indx++];
               const StorageSite::ScatterMap&  scatterMap = _meshList.at(id)->getCells().getScatterMap();
               const Array<int>& scatterArray = *(scatterMap.find(pos)->second);
               const int neighID = pos->getGatherProcID();
               _debugFile << " meshID = " << id << "   procID = " << _procID << "  neighProcID = " << neighID <<  " Tag = " << pos->getTag() << " : " << endl;
               for ( int i = 0; i < scatterArray.getLength(); i++ ){
                    _debugFile << "      scatterArray[" << i << "] = " << scatterArray[i]  << endl;
                }
           } else { //this means mesh interface
               const int id = scatterSiteVecID[indx++];
               const StorageSite::ScatterMap&  scatterMap = _meshList.at(id)->getCells().getScatterMap();
               const Array<int>& scatterArray = *(scatterMap.find(pos)->second);
               int neighMeshID = siteMeshMapper[pos];
               _debugFile << "   meshID = " << id << "   otherside MeshID = " << neighMeshID <<  " : " << endl;
               for ( int i = 0; i < scatterArray.getLength(); i++ ){
                   _debugFile << "      scatterArray[" << i << "] = " << scatterArray[i]  << endl;
               }
            }
     }
     debug_file_close();
}
void

MeshDismantler::debug_gather_mappers()
{

    debug_file_open("gatherMappers");
    //creating mappers between cell storage site to mesh id
    map< const StorageSite*, int > siteMeshMapper; //key  = storage site, value = mesh ID of cellSite
    for ( int id = 0; id < _nmesh; id++ )
        siteMeshMapper[&_meshList.at(id)->getCells()] = id;

    //first make sure consistent order for gatherMap since map order with key which is pointer
    //first mesh interface (increasing order) then partition interface (increasing order)
    vector< const StorageSite* > gatherSiteVec;
    map <int, int> packIDToMeshID;
    vector< int > gatherSiteVecID;
    IntStorageSiteMap   orderedMeshInterface;
    for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite::GatherMap&  gatherMap = _meshList.at(id)->getCells().getGatherMap();
        foreach ( const StorageSite::GatherMap::value_type& pos, gatherMap ){
            const StorageSite& gatherSite = *pos.first;
            if ( gatherSite.getGatherProcID() == -1 ){  //this means mesh interface
                int neighMeshID = siteMeshMapper[ pos.first ];
                const int packed_info = ( (id  << 16)  | neighMeshID );
                orderedMeshInterface.insert( make_pair<int, const StorageSite*>(packed_info, pos.first) );
                packIDToMeshID.insert( make_pair<int,int>(packed_info,id) );
            }
        }
    }

   //fill gatherSiteVec, map data structure already ordered them
   foreach ( const IntStorageSiteMap::value_type& pos, orderedMeshInterface ){
          gatherSiteVec.push_back( pos.second );
          const int meshID = packIDToMeshID[pos.first];
          gatherSiteVecID.push_back( meshID );
  }
   //clear for usage in partition interface
    packIDToMeshID.clear();
    //order partition interface (map will reorder)
    IntStorageSiteMap orderedPartInterface;
    for ( int id = 0; id < _nmesh; id++ ){
        const StorageSite::GatherMap&  gatherMap = _meshList.at(id)->getCells().getGatherMap();
        foreach ( const StorageSite::GatherMap::value_type& pos, gatherMap ){
            const StorageSite& gatherSite = *pos.first;
            if ( gatherSite.getGatherProcID() != -1 ){  //this means partition face
                const int neighID = gatherSite.getGatherProcID();
                const int tag     = gatherSite.getTag();
                const int packed_info = ( (id << 31)  | (neighID  << 24) | tag );
                orderedPartInterface.insert( make_pair<int, const StorageSite*>(packed_info, pos.first) );
                packIDToMeshID.insert( make_pair<int,int>(packed_info,id) );
            }
        }
     }

    //fill scatterSiteVec, map data structure already 
    foreach ( const IntStorageSiteMap::value_type& pos, orderedPartInterface ){
          gatherSiteVec.push_back( pos.second );
          const int meshID = packIDToMeshID[pos.first];
          gatherSiteVecID.push_back( meshID );
    }

     int indx = 0;
     foreach( const vector<const StorageSite*>::value_type& pos, gatherSiteVec){
          if ( pos->getGatherProcID() != -1 ){  //this means partition face
               const int id = gatherSiteVecID[indx++];
               const StorageSite::GatherMap&  gatherMap = _meshList.at(id)->getCells().getGatherMap();
               const Array<int>& gatherArray = *(gatherMap.find(pos)->second);
               const int neighID = pos->getGatherProcID();
               _debugFile << " meshID = " << id << "   procID = " << _procID << "  neighProcID = " << neighID <<  " Tag = " << pos->getTag() << " : " << endl;
               for ( int i = 0; i < gatherArray.getLength(); i++ ){
                    _debugFile << "      gatherArray[" << i << "] = " << gatherArray[i]  << endl;
                }
           } else { //this means mesh interface
               const int id = gatherSiteVecID[indx++];
               const StorageSite::GatherMap&  gatherMap = _meshList.at(id)->getCells().getGatherMap();
               const Array<int>& gatherArray = *(gatherMap.find(pos)->second);
               int neighMeshID = siteMeshMapper[pos];
               _debugFile << "   meshID = " << id << "   otherside MeshID = " << neighMeshID <<  " : " << endl;
               for ( int i = 0; i < gatherArray.getLength(); i++ ){
                   _debugFile << "      gatherArray[" << i << "] = " << gatherArray[i]  << endl;
               }
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
