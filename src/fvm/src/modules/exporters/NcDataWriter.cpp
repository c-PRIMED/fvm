//
// C++ Implementation: NcDataWriter
//
// Description: 
//
//
// Author: yildirim,,, <yildirim@prism>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <iostream>
#include <cassert>

#include "NcDataWriter.h"
#include "netcdfcpp.h"
#include "StorageSite.h"
#include "mpi.h"

NcDataWriter::NcDataWriter( const MeshList& meshes, const string& fname )
: _meshList( meshes ), _fname( fname ), MAX_CHAR(40)
{

   init();
   cout << " input file name = " << _fname << endl;
}


NcDataWriter::~NcDataWriter()
{
   if ( _ncFile ) delete _ncFile;
// // // //    vector< char* > ::iterator it;
// // // //    for ( it = _boundaryTypeVals.begin(); it != _boundaryTypeVals.end(); it++)
// // // //        delete [] *it;

}


void  
NcDataWriter::record()
{
    setNcFile();
    setDims();
    setVars();
   set_var_values();

   _ncFile->close();
}






                          //PRIVATE FUNCTIONS

void
NcDataWriter::init()
{
    _ncFile = NULL;

}

//setting NcFile
void
NcDataWriter::setNcFile()
{
    assert( !_ncFile );
    _ncFile = new NcFile( _fname.c_str(), NcFile::Replace );
    assert ( _ncFile->is_valid() );

}

//set NcDims
void
NcDataWriter::setDims()
{
    _nmesh = _ncFile->add_dim("nmesh",  long( _meshList.size() ) );
     
    int index = 0;
     for ( int id = 0; id < _nmesh->size(); id++ ){
         index += _meshList.at(id)->getBoundaryFaceGroups().size();
     }

      _nBoun  = _ncFile->add_dim("boun_type_dim", index );
      _charSize= _ncFile->add_dim("char_dim", MAX_CHAR );
      assert( _ncFile->add_att("nmesh", "number of total meshes") );
      assert( _ncFile->add_att("boun_type_dim", "total count of boundary types") );
      assert( _ncFile->add_att("char_dim", "maximum capacity of char variable" ) );

}

//set  NcVars
void 
NcDataWriter::setVars()
{

    _dimension =  _ncFile->add_var("dimension", ncInt, _nmesh );
    _meshID    =  _ncFile->add_var("mesh_id"  , ncInt, _nmesh );
    _facesCount = _ncFile->add_var("faces_count", ncInt, _nmesh );
    _cellsCount = _ncFile->add_var("cells_count", ncInt, _nmesh );
    _ghostCellsCount = _ncFile->add_var("ghost_cells_count", ncInt, _nmesh);
    _nodesCount = _ncFile->add_var("nodes_count", ncInt, _nmesh );
    _mapCount   = _ncFile->add_var("map_count", ncInt, _nmesh );
    _interiorFaceGroup = _ncFile->add_var("interior_faces_group", ncInt, _nmesh);
    _boundaryGroup  = _ncFile->add_var("boundary_group", ncInt, _nmesh); 
    _boundarySize   = _ncFile->add_var("boundary_size", ncInt, _nBoun);
    _boundaryOffset = _ncFile->add_var("boundary_offset", ncInt, _nBoun);
    _boundaryID     = _ncFile->add_var("boundary_id", ncInt, _nBoun);
    _boundaryType   = _ncFile->add_var("boundary_type", ncChar, _nBoun, _charSize );

}

//assign values to NcVars
void 
NcDataWriter::set_var_values()
{
    for ( long id = 0; id < _nmesh->size(); id++ ){
       //dimension-s
       _dimensionVals.push_back( _meshList.at(id)->getDimension() );

       //mesh id-s
       _meshIDVals.push_back( _meshList.at(id)->getID() );

       //faces, cell, and node counts
       _facesCountVals.push_back( _meshList.at(id)->getFaces().getSelfCount() );
       _cellsCountVals.push_back( _meshList.at(id)->getCells().getSelfCount() );
       _ghostCellsCountVals.push_back( _meshList.at(id)->getCells().getCount() - 
                                       _meshList.at(id)->getCells().getSelfCount() );
       _nodesCountVals.push_back( _meshList.at(id)->getNodes().getSelfCount() );

      //neighbour counts
       const Mesh::GhostCellSiteMap&  ghostCellSiteMap = _meshList.at(id)->getGhostCellSiteMap();
       _mapCountVals.push_back(  ghostCellSiteMap.size() );

        //interior face counts
       _interiorFaceGroupVals.push_back( _meshList.at(id)->getInteriorFaceGroup().site.getCount() );

       //boundary face
       const FaceGroupList& bounFaceList = _meshList.at(id)->getBoundaryFaceGroups();
       _boundaryGroupVals.push_back ( bounFaceList.size() );
       for ( int boun = 0; boun < int(bounFaceList.size()); boun++ ) {
           _boundarySizeVals.push_back( bounFaceList.at(boun)->site.getCount() );
           _boundaryOffsetVals.push_back( bounFaceList.at(boun)->site.getOffset() );
           _boundaryIDVals.push_back( bounFaceList.at(boun)->id );
           _boundaryTypeVals.push_back( bounFaceList.at(boun)->groupType.c_str() );

           _boundaryType->set_cur(boun); 
           assert(  int(bounFaceList.at(boun)->groupType.size()) < MAX_CHAR );
           _boundaryType->put( _boundaryTypeVals[boun], 1, bounFaceList.at(boun)->groupType.size() );
                      
          
       }



    }

     //add attributies
     assert( _dimension->add_att("dim", "dimension of meshes, 1:1D, 2:2D, 3:3D." ) );
     assert( _meshID->add_att("id", " mesh identificaton index" ) );
     assert( _facesCount->add_att("StorageSite", "number of faces") );
     assert( _cellsCount->add_att("StorageSite", "number of cells ") );
     assert( _ghostCellsCount->add_att("StorageSite", "number of ghost cells") );
     assert( _nodesCount->add_att("StorageSite", "number of nodes") );
     assert( _mapCount->add_att("neigh_count", "total neighboorhood mesh counts") );
     assert( _interiorFaceGroup->add_att("interior_face_group", "total interior faces") );
     assert( _boundaryGroup->add_att("boundary_group", " total boundary faces") );
     assert( _boundarySize->add_att("boundary_size", " size of boundary" ) );
     assert( _boundaryOffset->add_att("boundary_offset", " offset of boundary" ) );
     assert( _boundaryID->add_att("boundary_id", " boundary id " ) );
     assert( _boundaryType->add_att("boundary_type", " type of boundary condition ") );


     //write values
    _dimension->put( &_dimensionVals[0], _nmesh->size()  );
    _meshID->put( &_meshIDVals[0], _nmesh->size() );
    _facesCount->put( &_facesCountVals[0], _nmesh->size() );
    _cellsCount->put( &_cellsCountVals[0], _nmesh->size() );
    _ghostCellsCount->put( &_ghostCellsCountVals[0], _nmesh->size() );
    _nodesCount->put( &_nodesCountVals[0], _nmesh->size() );
    _mapCount->put( &_mapCountVals[0], _nmesh->size() );
    _interiorFaceGroup->put( &_interiorFaceGroupVals[0], _nmesh->size() );
    _boundaryGroup->put(&_boundaryGroupVals[0], _nmesh->size() );

    _boundarySize->put( &_boundarySizeVals[0], _boundarySizeVals.size()  );
    _boundaryOffset->put( &_boundaryOffsetVals[0], _boundaryOffsetVals.size() );
    _boundaryID->put( &_boundaryIDVals[0], _boundaryIDVals.size() );
    
      
    


}
