// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "MeshAssembler.h"
#include "CRConnectivity.h"
#include "MultiField.h"

using namespace std;

MeshAssembler::   MeshAssembler( const MeshList& meshList ):_meshList( meshList )
{

    init();
}

MeshAssembler::~MeshAssembler()
{
   vector< Mesh* >::iterator it_mesh;
   for ( it_mesh = _mesh.begin(); it_mesh != _mesh.end(); it_mesh++)
        delete *it_mesh;
}


void
MeshAssembler::init()
{
  int dim = _meshList.at(0)->getDimension();
   //construct merged Linearsystem
  _mesh.push_back( new Mesh( dim) );


  _interfaceNodesSet.resize( _meshList.size() );
  _localInterfaceNodesToGlobalMap.resize( _meshList.size() );

   setCellsSite();
   setFacesSite();
   setInterfaceNodes();
   //countInterfaceNodes();
   setNodesSite();
   setCellsMapper();
   setFaceCells();
   setNodesMapper();
   setFaceNodes();
   setCoord();
   setMesh();
   setMeshCellColor();

}

//setting storagesite for cells
void 
MeshAssembler::setCellsSite()
{
   int siteCount     = 0;
   int siteSelfCount = 0;
   for ( unsigned int id = 0; id < _meshList.size(); id++ ){
      const StorageSite& site = _meshList[id]->getCells();
      const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      int nghost = 0;
      foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
         const Array<int>& fromIndices = *(mpos.second);
         nghost += fromIndices.getLength();
      }
      siteCount += site.getCount() - nghost;
      siteSelfCount += site.getSelfCount();
   }

   _cellSite = StorageSitePtr( new StorageSite(siteSelfCount, siteCount-siteSelfCount) );

}

//setting storagesite for faces
void 
MeshAssembler::setFacesSite()
{
   int faceCount = 0;
   for ( unsigned int id = 0; id < _meshList.size(); id++ ){
      const StorageSite& site = _meshList[id]->getFaces();
      faceCount += site.getCount();
   }

   int sharedFaceCount = 0;
   for ( unsigned int id = 0; id < _meshList.size(); id++ ){
        const FaceGroupList &  faceGroupList = _meshList[id]->getInterfaceGroups();
       for ( int n = 0; n < _meshList[id]->getInterfaceGroupCount(); n++ ){
          const StorageSite& site = faceGroupList[n]->site;
          sharedFaceCount += site.getCount();
       }
   }
   _faceSite = StorageSitePtr( new StorageSite(faceCount - sharedFaceCount / 2 ) );
    assert( sharedFaceCount%2  == 0 );

}


void
MeshAssembler::setInterfaceNodes()
{
   //loop over meshes
   for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      const Mesh& mesh = *(_meshList.at(n));
      const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
      //const CRConnectivity& faceCells = mesh.getAllFaceCells();	
      const FaceGroupList& faceGroupList = mesh.getInterfaceGroups();
      //looop over  interfaces
      for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
           int id = faceGroupList[i]->id;
           set<int>&  nodes   = _interfaceNodesSet[n][id];
           //loop over faces to fill nodes
           const StorageSite& site = faceGroupList[i]->site;
           int jstart =   site.getOffset(); 
           int jend   = jstart + site.getCount();
           for ( int face = jstart; face < jend; face++ ){
               //loop over face surronding nodes
               for ( int node = 0; node < faceNodes.getCount(face); node++ ){
                   nodes.insert( faceNodes(face,node) );
/*                   cout << " nemsh = " << n << "faceNodes(" << face << "," << node << ") = "  << faceNodes(face,node) << 
                           " faceCells(" << face << "," << "0) = " << faceCells(face,0) << 
                           " faceCells(" << face << "," << "1) = " << faceCells(face,1) << endl;*/
               }
           }
      }
   }

}

//this method counts interface nodes after merging, 
//it seems only way to match is using its coordinates 
void 
MeshAssembler::countInterfaceNodes()
{
   
  for ( unsigned int n = 0; n < _meshList.size(); n++ ){
     const Mesh& mesh = *(_meshList.at(n));
     const StorageSite& faceSite = mesh.getFaces();
     const CRConnectivity& cellFaces = mesh.getCellFaces();
     for ( int i = 0; i < faceSite.getCount(); i++ ){
         cout << " cellFaces[" << i << " ] = ";
         for ( int j = 0; j < cellFaces.getCount(i); j++ )
             cout << "  " << cellFaces(i,j);
         cout << endl;
     }
  }
  //writing mapping
  for ( unsigned int i = 0; i < _meshList.size(); i++ ){
       const StorageSite& thisSite = _meshList[i]->getCells();
       const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMap();
       foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
          const StorageSite& oSite    = *(mpos.first);
          const Array<int>& fromIndices = *(mpos.second);
          const Array<int>& toIndices   = *(oSite.getGatherMap().find(&thisSite)->second);
          cout << " fromIndices.getLength() = " << fromIndices.getLength() << endl;
          cout << " scatterArray " << endl;
          for ( int n  = 0; n < fromIndices.getLength(); n++){
             cout << "      fromIndices[" << n << "] = " << fromIndices[n] << "      toIndices[" << n << "] = " << toIndices[n] << endl;
          }
       }
   }

}

//setting Storage site for nodes
void 
MeshAssembler::setNodesSite()
{
  //count inner nodes for assembly
  int nodeCount       = getInnerNodesCount();
  int nInterfaceNodes = getInterfaceNodesCount();
  _nodeSite = StorageSitePtr( new StorageSite(nodeCount + nInterfaceNodes) );

}

//gettin localToGlobal and globalToLocal for cell
void
MeshAssembler::setCellsMapper()
{
    //count cells
    int selfCount = 0;
    for ( unsigned int id = 0; id < _meshList.size(); id++ ){
       const StorageSite& cellSite = _meshList[id]->getCells();
       selfCount += cellSite.getSelfCount();
    }
    _globalCellToMeshID.resize( selfCount );
    _globalCellToLocal.resize ( selfCount );
    //loop over meshes
    int glblIndx = 0;
    for ( unsigned int id = 0; id < _meshList.size(); id++ ){
       const StorageSite& cellSite = _meshList[id]->getCells();
       _localCellToGlobal[id]   = ArrayIntPtr( new Array<int>( cellSite.getCount() ) ); 
        Array<int>&  localToGlobal = *_localCellToGlobal[id];
        localToGlobal = -1; //initializer 
       //loop over inner cells
       for ( int n = 0; n < cellSite.getSelfCount(); n++ ){
            localToGlobal[n] = glblIndx;
           _globalCellToMeshID[glblIndx] = id;  //belongs to which mesh from global Cell ID
           _globalCellToLocal [glblIndx] = n;   //gives local Cell ID from Global ID but which mesh comes from mapGlobalCellToMeshID
            glblIndx++;
       }
    }

  //creating cellID MultiField to use sync() operation
  shared_ptr<MultiField>  cellMultiField = shared_ptr<MultiField>( new MultiField()    );
  shared_ptr<Field>       cellField      = shared_ptr<Field>     ( new Field("cellID") );
 
  for ( unsigned int n = 0; n < _meshList.size(); n++ ){
     const StorageSite* site = &_meshList[n]->getCells();
     MultiField::ArrayIndex ai( cellField.get(), site );
     shared_ptr<Array<int> > cIndex(new Array<int>(site->getCount()));
      *cIndex = -1;
     cellMultiField->addArray(ai,cIndex);
  }

  //fill  local mesh with global Indx
  for ( unsigned int n = 0; n < _meshList.size(); n++ ){
     const StorageSite* site = &_meshList[n]->getCells();
     MultiField::ArrayIndex ai( cellField.get(), site );
     Array<int>&  localCell = dynamic_cast< Array<int>& >( (*cellMultiField)[ai] ); 
     const Array<int>&  localToGlobal = *_localCellToGlobal[n];
     for ( int i = 0; i < site->getSelfCount(); i++ )
          localCell[i] =localToGlobal[i];
  }
  //sync opeartion
  cellMultiField->sync();

  //fill  local mesh with global Indx after sync
  for ( unsigned int n = 0; n < _meshList.size(); n++ ){
     const StorageSite* site = &_meshList[n]->getCells();
     MultiField::ArrayIndex ai( cellField.get(), site );
     const Array<int>&  localCell = dynamic_cast< const Array<int>& >( (*cellMultiField)[ai] ); 
     Array<int>&  localToGlobal = *_localCellToGlobal[n];
     for ( int i = 0; i < site->getCount(); i++ )
          localToGlobal[i] = localCell[i];
  }


//    //above algorithm fill localCellToGlobal for only inner cells but we will do ghost cells for interfaces
//    for ( unsigned int i = 0; i < _meshList.size(); i++ ){
//        const StorageSite& thisSite = _meshList[i]->getCells();
//        const StorageSite::ScatterMap& scatterMap = thisSite.getScatterMap();
//        foreach(const StorageSite::ScatterMap::value_type& mpos, scatterMap){
//           const StorageSite& oSite    = *(mpos.first);
//           const Array<int>& fromIndices = *(mpos.second);
//           const Array<int>& toIndices   = *(oSite.getGatherMap().find(&thisSite)->second);
//           cout << " fromIndices.getLength() = " << fromIndices.getLength() << endl;
//           cout << " scatterArray " << endl;
//           for ( int n  = 0; n < fromIndices.getLength(); n++){
//              cout << "      fromIndices[" << n << "] = " << fromIndices[n] << "      toIndices[" << n << "] = " << toIndices[n] << endl;
//           }
//        }
//    }

}

//getting CRConnectivity faceCells
void
MeshAssembler::setFaceCells()
{
    _faceCells = CRConnectivityPtr( new CRConnectivity( *_faceSite, *_cellSite) );
    _faceCells->initCount();
    
     //addCount
     const int cellCount = 2;
     for ( int i = 0; i < _faceSite->getCount(); i++ )
        _faceCells->addCount(i, cellCount); // face has always two cells
     //finish count
     _faceCells->finishCount();

     //first interior faces 
     int face = 0;
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
        const Mesh& mesh = *(_meshList[n]);
        const CRConnectivity& faceCells = mesh.getAllFaceCells();
        const FaceGroup& faceGroup = mesh.getInteriorFaceGroup();
        const Array<int>& localToGlobal = *_localCellToGlobal[n];
        for ( int i = 0; i < faceGroup.site.getCount(); i++ ){
            int cell1 = faceCells(i,0);
            int cell2 = faceCells(i,1);
            _faceCells->add( face, localToGlobal[ cell1 ] );
            _faceCells->add( face, localToGlobal[ cell2 ] );
            face++;
        }
     }

     //now add interfaces
     set<int> faceGroupSet; // this set make sure that no duplication in sweep
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
          const Mesh& mesh = *(_meshList[n]);
          const CRConnectivity& faceCells = mesh.getAllFaceCells();
          const FaceGroupList&  faceGroupList = mesh.getInterfaceGroups();
          const Array<int>& localToGlobal = *_localCellToGlobal[n];
          //loop over interfaces 
          for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
               int faceGroupID = faceGroupList[i]->id;
               //pass only if it is not in faceGroupSet
               if ( faceGroupSet.count( faceGroupID) == 0 ){
                   faceGroupSet.insert( faceGroupID );
                   //sweep interfaces for add up operation
                   int ibeg = faceGroupList[i]->site.getOffset();
                   int iend = ibeg + faceGroupList[i]->site.getCount();
                   for ( int i = ibeg; i < iend; i++ ){
                      int cell1 = faceCells(i,0);
                      int cell2 = faceCells(i,1);
                      _faceCells->add( face, localToGlobal[ cell1 ] );
                      _faceCells->add( face, localToGlobal[ cell2 ] );
                      face++;
                   }
               }
           }
      }
     //now add boundary faces
     int indx = _cellSite->getSelfCount();
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
          const Mesh& mesh = *(_meshList[n]);
          const StorageSite& cellSite = mesh.getCells();
          const CRConnectivity& faceCells = mesh.getAllFaceCells();
          const FaceGroupList&  bounGroupList = mesh.getBoundaryFaceGroups();
          const Array<int>& localToGlobal = *_localCellToGlobal[n];
          //loop over interfaces 
          for ( int i = 0; i <mesh.getBoundaryGroupCount(); i++ ){
              //sweep interfaces for add up operation
              int ibeg = bounGroupList[i]->site.getOffset();
              int iend = ibeg + bounGroupList[i]->site.getCount();
              for ( int i = ibeg; i < iend; i++ ){
                  int cell1 = faceCells(i,0);
                  int cell2 = faceCells(i,1);
                  //cell1
                  if ( cell1 < cellSite.getSelfCount() )
                      _faceCells->add( face, localToGlobal[ cell1 ] ); 
                  else 
                      _faceCells->add( face, localToGlobal[ cell2 ] ); 
                  //adding boundary cells 	
                  _faceCells->add( face, indx  ); 
                  indx++;
                  face++;
              }
          }
       }

      //finishAdd
      _faceCells->finishAdd();


}

//getting CRConnectivity faceNodes
void
MeshAssembler::setFaceNodes()
{
    _faceNodes = CRConnectivityPtr( new CRConnectivity( *_faceSite, *_nodeSite) );
    _faceNodes->initCount();
     //addCount
     const Mesh& mesh0 = *_meshList[0];
     const int nodeCount = mesh0.getAllFaceNodes().getRow()[1] -  mesh0.getAllFaceNodes().getRow()[0];
     for ( int i = 0; i < _faceSite->getCount(); i++ )
        _faceNodes->addCount(i, nodeCount); 
     //finish count
     _faceNodes->finishCount();
     //first interior faces 
     int face = 0;
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
        const Mesh& mesh = *(_meshList[n]);
        const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
        const FaceGroup& faceGroup = mesh.getInteriorFaceGroup();
        const StorageSite& faceGroupSite = faceGroup.site;
        const Array<int>& localToGlobal = *_localNodeToGlobal[n];
        for ( int i = 0; i < faceGroupSite.getCount(); i++ ){
           for ( int j = 0; j < faceNodes.getCount(i); j++ ){
               int nodeID = faceNodes(i,j);
              _faceNodes->add( face, localToGlobal[nodeID] );
           }
           face++;
        }
     }

     //interfaces (mesh)
     set<int> faceGroupSet;
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
         const Mesh& mesh = *(_meshList[n]);
         const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
         const FaceGroupList& faceGroupList = mesh.getInterfaceGroups();
         const Array<int>& localToGlobal = *_localNodeToGlobal[n];
         for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
               int faceGroupID = faceGroupList[i]->id;
               //pass only if it is not in faceGroupSet
               if ( faceGroupSet.count( faceGroupID) == 0 ){
                   faceGroupSet.insert( faceGroupID );
                   //sweep interfaces for add up operation
                   int ibeg = faceGroupList[i]->site.getOffset();
                   int iend = ibeg + faceGroupList[i]->site.getCount();
                   for ( int i = ibeg; i < iend; i++ ){
                       for ( int j = 0; j < faceNodes.getCount(i); j++ ){
                           int nodeID = faceNodes(i,j);
                           _faceNodes->add( face, localToGlobal[nodeID] );
                       }
                       face++;
                   }
                }
         }
     }
     //interior face size stored
     _interiorFaceSize = face;


     //now add boundary faces
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
          const Mesh& mesh = *(_meshList[n]);
          const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
          const FaceGroupList&  bounGroupList = mesh.getBoundaryFaceGroups();
          const Array<int>& localToGlobal = *_localNodeToGlobal[n];
          //loop over interfaces 
          for ( int i = 0; i <mesh.getBoundaryGroupCount(); i++ ){
              //sweep interfaces for add up operation
              int ibeg = bounGroupList[i]->site.getOffset();
              int iend = ibeg + bounGroupList[i]->site.getCount();
              for ( int i = ibeg; i < iend; i++ ){
                  for ( int j = 0; j < faceNodes.getCount(i); j++ ){
                      int nodeID = faceNodes(i,j);
                      _faceNodes->add( face, localToGlobal[nodeID] );
                   }
                  face++;
              }
          }
       }

    //finishAdd
    _faceNodes->finishAdd();
}

//count only inner nodes of all mesh after assembly
int
MeshAssembler::getInnerNodesCount(){  
   int nodeCount = 0;
   for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      const Mesh& mesh = *(_meshList.at(n));
      const CRConnectivity& faceNodes = mesh.getAllFaceNodes();
      const StorageSite& site = _meshList[n]->getNodes();
      Array<int>  nodeArray( site.getCount() );
      nodeArray = -1;
      //looop over  interior faces
      int nface = mesh.getFaces().getCount();
      for ( int i = 0; i < nface; i++ )
          //loop over face surronding nodes
           for ( int node = 0; node < faceNodes.getCount(i); node++ )
               nodeArray[ faceNodes(i,node) ] = 1;

    //loop over interface faces to reset 
     const FaceGroupList& faceGroupList = mesh.getInterfaceGroups();  
     for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
           int id = faceGroupList[i]->id;
           set<int>&  nodes   = _interfaceNodesSet[n][id];
           foreach ( const set<int>::value_type node, nodes )
               nodeArray[ node] = -1;  //reseting again, we want to seperate inner/boundary nodes than interface nodes
     }
     //node count
     for ( int i = 0; i < nodeArray.getLength(); i++ )
         if ( nodeArray[i] != -1 ) nodeCount++;
   }

  return nodeCount;
}

//count nodes (duplicated) on shared interfaces from each local mesh side
int
MeshAssembler::getInterfaceNodesDuplicatedCount()
{
   //loop over meshes
   int nInterfaceNodes = 0;
   for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      const Mesh& mesh = *(_meshList.at(n));
      const FaceGroupList& faceGroupList = mesh.getInterfaceGroups();
      //looop over  interfaces
      for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
           int id = faceGroupList[i]->id;
           const set<int>&  nodes   = _interfaceNodesSet[n][id];
           nInterfaceNodes += nodes.size();
      }
   }
   return nInterfaceNodes;
}

//count nodes on interfaces (not duplicated) 
int
MeshAssembler::getInterfaceNodesCount()
{
  //count nodes on shared interfaces (duplicated)
  int nInterfaceNodes = getInterfaceNodesDuplicatedCount();

   Array< Mesh::VecD3 > interfaceNodeValues ( nInterfaceNodes );
   Array< int > globalIndx  ( nInterfaceNodes );
   globalIndx   = -1;

   //filing interfaceNodeValues
   int indx = 0;
   for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      const Mesh& mesh = *(_meshList.at(n));
      const FaceGroupList& faceGroupList = mesh.getInterfaceGroups();  
      const Array<Mesh::VecD3>& coord = mesh.getNodeCoordinates();
      //looop over  interfaces
      for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
           int id = faceGroupList[i]->id;
           const set<int>&  nodes   = _interfaceNodesSet[n][id];
           foreach ( const set<int>::value_type node, nodes ){
              interfaceNodeValues[indx] = coord[node];
              indx++;
           }
      }
   }

   //greedy algorithm to fill interfaceNodeCount
  indx =0;
  for ( int i = 0; i < nInterfaceNodes; i++ ){
     if ( globalIndx[i] == -1 ){
        globalIndx[i] = indx;
        double x = interfaceNodeValues[i][0];
        double y = interfaceNodeValues[i][1];
        double z = interfaceNodeValues[i][2];
        for ( int j = i+1; j < nInterfaceNodes; j++ ){
           if ( globalIndx[j] == -1 ){
              double xOther = interfaceNodeValues[j][0];
              double yOther = interfaceNodeValues[j][1];
              double zOther = interfaceNodeValues[j][2];
              if ( x == xOther && y == yOther && z == zOther )
                 globalIndx[j] = indx;
           }
        }
        indx++;
     }
  }


  _nInterfaceNodes = indx;
  //filling localInterfaceNodes to GlobalNodes data structure
   indx = 0;
   for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      const Mesh& mesh = *(_meshList.at(n));
      const FaceGroupList& faceGroupList = mesh.getInterfaceGroups();
      //looop over  interfaces
      for ( int i = 0; i < mesh.getInterfaceGroupCount(); i++ ){
          int id = faceGroupList[i]->id;
          const set<int>&  nodes   = _interfaceNodesSet[n][id];
          foreach ( const set<int>::value_type node, nodes ){
             _localInterfaceNodesToGlobalMap[n][node] = globalIndx[indx];
             indx++;
          }
      }
   }

   return _nInterfaceNodes;
}

//local Nodes to global Nodes 
void
MeshAssembler::setNodesMapper()
{

    int glblIndx = _nInterfaceNodes; // node numbering first from interfaces
    cout << " glbIndx = " << glblIndx << endl;
    for ( unsigned int id = 0; id < _meshList.size(); id++ ){
       const StorageSite& nodeSite = _meshList[id]->getNodes();
       const StorageSite& faceSite = _meshList[id]->getFaces();
       const CRConnectivity& faceNodes = _meshList[id]->getAllFaceNodes();
       _localNodeToGlobal[id]   = ArrayIntPtr( new Array<int>( nodeSite.getCount() ) ); 
       Array<bool> isVisitedNodes( nodeSite.getCount() );
       isVisitedNodes = false;
        Array<int>&  localToGlobal = *_localNodeToGlobal[id];
        localToGlobal = -1; //initializer 
        const  map<int,int>&  localInterfaceNodesToGlobalMap = _localInterfaceNodesToGlobalMap[id];
        //loop over faces then connecting nodes
        for ( int face = 0; face < faceSite.getCount(); face++ ){
           for ( int n = 0; n < faceNodes.getCount(face); n++ ){
               int nodeID = faceNodes(face,n); 
               //cout << " nodeID = " << nodeID << " isVisitedNodes = " << isVisitedNodes[nodeID] << endl;
               if ( !isVisitedNodes[nodeID]){ //if it is not visited, we renumber it
                   if ( localInterfaceNodesToGlobalMap.count( nodeID ) == 0 ) 
                      localToGlobal[nodeID] = glblIndx++;
                   else 
                      localToGlobal[nodeID] = localInterfaceNodesToGlobalMap.find(nodeID)->second;
                   isVisitedNodes[nodeID] = true;
               }
           }
        }
    }



}


void
MeshAssembler::setBoundaryFaceGroup()
{
     //interior face size stored
     int face =   _interiorFaceSize;
     //now add boundary faces
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
          const Mesh& mesh = *(_meshList[n]);
          const FaceGroupList&  bounGroupList = mesh.getBoundaryFaceGroups();
          //loop over interfaces 
          for ( int i = 0; i <mesh.getBoundaryGroupCount(); i++ ){
              //sweep interfaces for add up operation
              int size   = bounGroupList[i]->site.getCount();
              int offset = face;
              int id     = bounGroupList[i]->id;
              string  boundaryType = bounGroupList[i]->groupType;
              //update face for offset
              face += size;
              _mesh.at(0)->createBoundaryFaceGroup( size, offset, id, boundaryType );
          }
     }

}

//setting coordinates
void
MeshAssembler::setCoord()
{
   _coord = ArrayVecD3Ptr(new Array<Mesh::VecD3>(_nodeSite->getCount()));
    for ( unsigned n = 0; n < _meshList.size(); n++ ){
         const Mesh& mesh = *(_meshList[n]);
         const StorageSite& nodeSite = mesh.getNodes();
         const Array<int>& localToGlobal = *_localNodeToGlobal[n];
         const Array<Mesh::VecD3>&   coordLocal = mesh.getNodeCoordinates();
         //loop over local mesh nodes
         for ( int i = 0; i < nodeSite.getCount(); i++ ){
               if ( localToGlobal[i] != -1 ) //so mapping is well defined
                 (*_coord)[ localToGlobal[i] ] = coordLocal[i];
         }
    }
}

//set StorageSites
void
MeshAssembler::setSites()
{
   StorageSite& faceSite = _mesh.at(0)->getFaces();
   StorageSite& cellSite = _mesh.at(0)->getCells();
   StorageSite& nodeSite = _mesh.at(0)->getNodes();
   faceSite.setCount( _faceSite->getCount() );
   cellSite.setCount( _cellSite->getSelfCount(), _cellSite->getCount()-_cellSite->getSelfCount() );
   nodeSite.setCount( _nodeSite->getCount() );
}
//setting assembled mesh
void
MeshAssembler::setMesh()
{
    //set sites
    setSites();
    //interior face group
    _mesh.at(0)->createInteriorFaceGroup( _interiorFaceSize );
    //boundary face group
    setBoundaryFaceGroup();
    //setting coordinates
    _mesh.at(0)->setCoordinates( _coord     );
    //setting faceNodes CRConnecitivty
    _mesh.at(0)->setFaceNodes  ( _faceNodes );
    //setting faceCells CRConnectivity
    _mesh.at(0)->setFaceCells  ( _faceCells );

}

//filling _cellColor array in the merged mesh. This will be used in Parmetis (elmWghts)
void
MeshAssembler::setMeshCellColor()
{
     //allocate Mesh::_cellColor
    _mesh.at(0)->createCellColor();

    Array<int>& cellColor = _mesh.at(0)->getCellColors();
    //first interior faces of local meshes
    for ( unsigned int n = 0; n < _meshList.size(); n++ ){
       const Mesh& mesh = *(_meshList[n]);
       const CRConnectivity& faceCells = mesh.getAllFaceCells();
       const FaceGroup& faceGroup = mesh.getInteriorFaceGroup();
       const Array<int>& localToGlobal = *_localCellToGlobal[n];
       for ( int i = 0; i < faceGroup.site.getCount(); i++ ){
           int cell1 = faceCells(i,0);
           int cell2 = faceCells(i,1);
           cellColor[ localToGlobal[ cell1 ] ] = n;
           cellColor[ localToGlobal[ cell2 ] ] = n;
       }
    }

  //now color only boundary cells
     int indx = _cellSite->getSelfCount();
     for ( unsigned int n = 0; n < _meshList.size(); n++ ){
          const Mesh& mesh = *(_meshList[n]);
          const FaceGroupList&  bounGroupList = mesh.getBoundaryFaceGroups();
          //loop over interfaces 
          for ( int i = 0; i <mesh.getBoundaryGroupCount(); i++ ){
              //sweep interfaces for add up operation
              int ibeg = bounGroupList[i]->site.getOffset();
              int iend = ibeg + bounGroupList[i]->site.getCount();
              for ( int i = ibeg; i < iend; i++ ){
                  cellColor[indx] = n;
                  //adding boundary cells 	
                  indx++;
              }
          }
       }
      //assigning Mesh::numOfAssembleMesh
     _mesh[0]->setNumOfAssembleMesh( int(_meshList.size()) );

}


void 
MeshAssembler::debug_sites()
{
  debug_file_open("sites");
  _debugFile << " cells.getSelfCount() = " << _cellSite->getSelfCount() << " cells.selfCount() = " << _cellSite->getCount() << endl;
  _debugFile << " faces.getSelfCount() = " << _faceSite->getSelfCount() << " faces.selfCount() = " << _faceSite->getCount() << endl;
  _debugFile << " nodes.getSelfCount() = " << _nodeSite->getSelfCount() << " nodes.selfCount() = " << _nodeSite->getCount() << endl;
  debug_file_close();
}


void
MeshAssembler::debug_localToGlobal_mappers()
{
   debug_file_open("localToGlobal");
  //mapLocalToGlobal
   for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      _debugFile << " mesh = " << n << endl;
      const Array<int>& localToGlobal = *_localCellToGlobal[n];
      for ( int i = 0; i < localToGlobal.getLength(); i++ )
          _debugFile  << " localCellToGlobal[" << i << "] = " << localToGlobal[i] << endl;
      _debugFile << endl;
   }
   debug_file_close();
}

void 
MeshAssembler::debug_globalCellToMeshID_mappers()
{
   debug_file_open("globalCellToMeshID");
   //globalCellToMeshID
   for (  unsigned int i = 0; i < _globalCellToMeshID.size(); i++ )
      _debugFile  << " globalCellToMeshID[" << i << "] = " << _globalCellToMeshID[i] << endl;
  
   _debugFile << endl;
   //globalCellToLocal
   for ( unsigned int i = 0; i < _globalCellToLocal.size(); i++ )
      _debugFile  << " globalCellToLocal[" << i << "] = " << _globalCellToLocal[i] << endl;
    
   debug_file_close();
}


void
MeshAssembler::debug_sync_localToGlobal_mappers()
{
  debug_file_open("syncLocalToGlobal");
  //printing things
  _debugFile << " localCellToGlobal after sync() opeartion " << endl;
  for ( unsigned int n = 0; n < _meshList.size(); n++ ){
     const Array<int>&  localToGlobal = *_localCellToGlobal[n];
     _debugFile << " mesh = " << n << endl;
     for ( int i = 0; i < localToGlobal.getLength(); i++ )
         _debugFile << " localToGlobal[" << i << "] = " << localToGlobal[i] << endl;
     _debugFile << endl;
  }
  debug_file_close();
}


void
MeshAssembler::debug_faceCells()
{
  debug_file_open("faceCells");
  //faceCells
  _debugFile << " faceCells Connectivity " << endl;
  for ( int i = 0; i < _faceSite->getCount(); i++ ) { 
      int ncells = _faceCells->getCount(i);
      for ( int j = 0; j < ncells; j++ ){
         _debugFile << " faceCells(" << i << "," << j << ") = " << (*_faceCells)(i,j);
      }
      _debugFile << endl;
  }
  debug_file_close();
}

void
MeshAssembler::debug_localNodeToGlobal()
{
  debug_file_open("localNodeToGlobal");
  //localNodeToGlobal
  _debugFile << " localNodeToGlobal " << endl;
  for ( unsigned int n = 0; n < _meshList.size(); n++ ){
      const Array<int>&  localToGlobal = *_localNodeToGlobal[n];
      for ( int i = 0; i < localToGlobal.getLength(); i++ )
         _debugFile << " localToGlobal[" << i << "] = " << localToGlobal[i] << endl;
      _debugFile << endl;
  }
  debug_file_close();
}

void  
MeshAssembler::debug_print()
{
    debug_sites();
    debug_localToGlobal_mappers();
    debug_globalCellToMeshID_mappers();
    debug_sync_localToGlobal_mappers();
    debug_faceCells();
    debug_localNodeToGlobal();

}



void
MeshAssembler::debug_file_open( const string& fname_ )
{  
     string  fname = "MESHASSEMBLER_"+fname_+".dat";
    _debugFile.open( fname.c_str() );
}

void
MeshAssembler::debug_file_close()
{
    _debugFile.close();
}
