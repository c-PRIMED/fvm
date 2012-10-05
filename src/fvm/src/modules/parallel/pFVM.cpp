// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <string>
#include <vector>

#include "MMReader.h"
#include "AMG.h"
#include "Array.h"

#include "PartMesh.h"

#include "FluentReader.h"
#include "Mesh.h"
#include "NcDataWriter.h"
#include "NcDataReader.h"

using namespace std;

void  debug_mesh(MeshList mesh_list);


int main(int argc, char *argv[])
{
/*  MMReader reader(argv[1], argv[2]);
  
  shared_ptr<LinearSystem> ls(reader.getLS());

  AMG solver;
  solver.solve(*ls); */


   MPI::Init(argc, argv);

   //PartMesh*  partMesh = new PartMesh( MPI::COMM_WORLD.Get_rank() );


   string file_name( argv[1] );
   FluentReader* fluent_reader = new FluentReader( file_name );
   cout << " arg 1= " << argv[1] << endl;
   fluent_reader->readMesh();

   MeshList mesh_list = fluent_reader->getMeshList();

   //debug_mesh( mesh_list );


  //mesh_list
   vector<int> npart( mesh_list.size(),  MPI::COMM_WORLD.Get_size() );
   vector<int> etype( mesh_list.size(), 1);

  //constructer to PartMesh
   PartMesh*  part_mesh = new PartMesh(mesh_list, npart, etype );

  //set properties of Partition
   part_mesh->setWeightType(0);
   part_mesh->setNumFlag(0);

  //actions
//   part_mesh->fiedler_order("permutation26.txt");

   part_mesh->partition();

   part_mesh->mesh();
   
   delete fluent_reader;

//    stringstream ss;
//    ss << "test_" << MPI::COMM_WORLD.Get_rank() << ".cdf";
//    NcDataWriter  nc_writer( part_mesh.meshList(), ss.str() );
//    nc_writer.record();
// 
//    NcDataReader  nc_reader( ss.str() );
//    stringstream ss_test;
//    ss_test << "test_test_" << MPI::COMM_WORLD.Get_rank() << ".cdf";
//    MeshList meshes = nc_reader.getMeshList();
//    
//    NcDataWriter  nc_writer_test( meshes, ss_test.str() );
//    nc_writer_test.record();
//    NcDataReader::destroyMeshList( meshes );

   delete part_mesh;
   MPI::Finalize();

}





void  debug_mesh( MeshList mesh_list )
{
//understand Mesh 
   Mesh*    mesh = mesh_list.at(0);

   cout << " mesh dimension = " << mesh->getDimension() << endl;
   cout << " mesh ID        = " << mesh->getID()        << endl;
   
   StorageSite&  faces = mesh->getFaces();
   StorageSite&  cells = mesh->getCells();
   StorageSite&  nodes = mesh->getNodes();

   cout << " Faces " << endl;
   cout << " StorageSite::getCount()     = " << faces.getCount() << endl;
   cout << " StorageSite::getSelfCount() = " << faces.getSelfCount() << endl;
   cout << " StorageSite::getOffset()    = " << faces.getOffset() << endl;

   cout << " Cells " << endl;
   cout << " StorageSite::getCount()     = " << cells.getCount() << endl;
   cout << " StorageSite::getSelfCount() = " << cells.getSelfCount() << endl;
   cout << " StorageSite::getOffset()    = " << cells.getOffset() << endl;

   cout << " Nodes " << endl;
   cout << " StorageSite::getCount()     = " << nodes.getCount() << endl;
   cout << " StorageSite::getSelfCount() = " << nodes.getSelfCount() << endl;
   cout << " StorageSite::getOffset()    = " << nodes.getOffset() << endl;

{
   const CRConnectivity&  cellNodes= mesh->getCellNodes();
   //const CRConnectivity&  cellNodes= mesh->getCellCells();   
   //const CRConnectivity&    cellNodes = mesh->getCellFaces();
   //const CRConnectivity&   cellNodes = mesh->getAllFaceNodes();
   const Array<int>&   rowCellNodes = cellNodes.getRow();
   const Array<int>&   colCellNodes = cellNodes.getCol();
   

   cout << " CRConnectivity::getCellNodes() " << endl;
   cout << " ::getRowDim() = " << cellNodes.getRowDim() << endl;
   cout << " ::getColDim() = " << cellNodes.getColDim() << endl;
   for ( int cell = 0; cell < cellNodes.getRowDim(); cell++){
      cout << " row[" << cell+1 <<"] = " << rowCellNodes[cell] << "    ";
      int nnodes = rowCellNodes[cell+1] - rowCellNodes[cell];
      for ( int node = 0; node < nnodes; node++){
         cout << colCellNodes[ rowCellNodes[cell] + node ]+1 << "    ";
      }
      cout << endl;
   }
}
   cout << " " << endl;
   cout << " Mesh::getFaceGroupCount()     = " << mesh->getFaceGroupCount() << endl;
   cout << " Mesh::getBoundaryGroupCount() = " << mesh->getBoundaryGroupCount() << endl;
   cout << " Mesh::getInterfaceGroupCount()= " << mesh->getInterfaceGroupCount() << endl;
   cout << endl;
   
   const FaceGroup& interiorFaceGroup = mesh->getInteriorFaceGroup();
   cout << " interiorFaceGroup id        = " << interiorFaceGroup.id << endl;
   cout << " interiorFaceGroup groupType = " << interiorFaceGroup.groupType << endl; 
   cout << " site.getCount()             = " << interiorFaceGroup.site.getCount() << endl;
   cout << endl;
   
   const FaceGroupList&  bounFaceGroup = mesh->getBoundaryFaceGroups();

  
   
   for ( int i = 0; i < int(bounFaceGroup.size()); i++){
      cout << " bounFaceGroup id        = " << bounFaceGroup.at(i)->id << endl;
      cout << " bounFaceGroup groupType = " << bounFaceGroup.at(i)->groupType << endl; 
      cout << " site.getCount()             = " << bounFaceGroup.at(i)->site.getCount() << endl;
      cout << endl;
   }


{
  //use of CRConnectivity::getSubset
  const CRConnectivity&  cellNodes= mesh->getCellNodes();
  StorageSite site( 3, 2);
  Array<int> indices(5);
  indices[0] = 35; indices[1] = 3; indices[2] = 2; indices[3] = 16; indices[4] = 29;
  shared_ptr<CRConnectivity> subCellNodes = cellNodes.getSubset( site, indices);
  const Array<int>&   rowCellNodes = subCellNodes->getRow();
  const Array<int>&   colCellNodes = subCellNodes->getCol();

  cout << " subCellNodes = " << endl;
  cout << " getRowDim()  = " << subCellNodes->getRowDim() << endl;
  cout << " getColDim()  = " << subCellNodes->getColDim() << endl;
  for ( int cell = 0; cell < subCellNodes->getRowDim(); cell++){
      cout << " row[" << cell+1 <<"] = " << rowCellNodes[cell] << "    ";
      int nnodes = rowCellNodes[cell+1] - rowCellNodes[cell];
      for ( int node = 0; node < nnodes; node++){
         cout << colCellNodes[ rowCellNodes[cell] + node ]+1 << "    ";
      }
      cout << endl;
  }

 cout << endl;
}

//use of CRConnectivity::getLocalizedSubset(..)
{
  const CRConnectivity&  cellNodes= mesh->getCellNodes();
  StorageSite rowsite( 3, 2);
  StorageSite colsite(13);
  Array<int> indices(5);
  indices[0] = 35; indices[1] = 3; indices[2] = 2; indices[3] = 16; indices[4] = 29;
  shared_ptr<CRConnectivity> subLocalCellNodes = cellNodes.getLocalizedSubset( rowsite, colsite, indices);
  const Array<int>&   rowCellNodes = subLocalCellNodes->getRow();
  const Array<int>&   colCellNodes = subLocalCellNodes->getCol();

  cout << " subCellNodes = " << endl;
  cout << " getRowDim()  = " << subLocalCellNodes->getRowDim() << endl;
  cout << " getColDim()  = " << subLocalCellNodes->getColDim() << endl;
  for ( int cell = 0; cell < subLocalCellNodes->getRowDim(); cell++){
      cout << " row[" << cell+1 <<"] = " << rowCellNodes[cell] << "    ";
      int nnodes = rowCellNodes[cell+1] - rowCellNodes[cell];
      for ( int node = 0; node < nnodes; node++){
         cout << colCellNodes[ rowCellNodes[cell] + node ]+1 << "    ";
      }
      cout << endl;
  }

 cout << endl;

}


  
   


}
