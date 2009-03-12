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
#include <numeric>

#include "NcDataWriter.h"
#include "netcdfcpp.h"
#include "StorageSite.h"
#include "Array.h"
#include "OneToOneIndexMap.h"

int NcDataWriter::_writeAction = 0;

NcDataWriter::NcDataWriter( const MeshList& meshes, const string& fname )
: _meshList( meshes ), _fname( fname ), _ncFile(NULL), _xVals(NULL), _yVals(NULL), _zVals(NULL),
_fromIndicesVals(NULL), _toIndicesVals(NULL), MAX_CHAR(40), BOUN_TYPE_DIM(false), NEIGH_MESH( false ), INTERFACE( false )
{

   init();
}


NcDataWriter::~NcDataWriter()
{
   if ( _xVals ) delete [] _xVals;
   if ( _yVals ) delete [] _yVals;
   if ( _zVals ) delete [] _zVals;

   if ( _fromIndicesVals  ) delete []  _fromIndicesVals;
   if ( _toIndicesVals    ) delete []  _toIndicesVals;

   if ( _ncFile ) delete _ncFile;

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
   _writeAction++;
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
//set NcDims
void
NcDataWriter::setDims()
{

    _nmesh = _ncFile->add_dim("nmesh",   _meshList.size()  );
    assert( _ncFile->add_att("nmesh", "number of total meshes") );
    int index_boun = 0;
    int index_interface = 0;
    int  nnodes = 0;
    int  nfaces = 0;
    int  ncells = 0;
    int  nface_row = 0;
    int  nfaceCells_col = 0;
    int  nfaceNodes_col = 0; 
    int  ninterface = 0;

    //summing up values from all meshes
     for ( int id = 0; id < _nmesh->size(); id++ ){
         index_boun      += _meshList.at(id)->getBoundaryGroupCount();
         index_interface += _meshList.at(id)->getInterfaceGroupCount();
         nnodes          += _meshList.at(id)->getNodes().getCount();
         nfaces          += _meshList.at(id)->getFaces().getCount();
         ncells          += _meshList.at(id)->getCells().getCount();
         nface_row       += _meshList.at(id)->getAllFaceCells().getRow().getLength();
         nfaceCells_col  += _meshList.at(id)->getAllFaceCells().getCol().getLength();
         nfaceNodes_col  += _meshList.at(id)->getAllFaceNodes().getCol().getLength();

         const StorageSite::MappersMap&   mappers = _meshList.at(id)->getCells().getMappers();
         //loop over neighbour mesh to count interfaces
         StorageSite::MappersMap::const_iterator it;
         for ( it = mappers.begin(); it != mappers.end(); it++ )
             ninterface += it->second->getFromIndices().getLength();  //or it->secon->toIndices.getLength()

     }


      if ( index_boun > 0 ){
          BOUN_TYPE_DIM = true;
         _nBoun     = _ncFile->add_dim("boun_type_dim", index_boun );
         assert( _ncFile->add_att("boun_type_dim", "total count of boundary types") );
      }


      _charSize  = _ncFile->add_dim("char_dim", MAX_CHAR );
      assert( _ncFile->add_att("char_dim", "maximum capacity of char variable" ) );

      if ( index_interface > 0 ){
         NEIGH_MESH = true;
         _nNeighMesh= _ncFile->add_dim("nNeighMesh", index_interface ); 
          assert( _ncFile->add_att("nNeighMesh", "count of neighbour meshes")  );
      }

      _nnodes     = _ncFile->add_dim("nnodes", nnodes);
      _nfaces     = _ncFile->add_dim("nfaces", nfaces);
      _ncells     = _ncFile->add_dim("ncells", ncells);
      assert( _ncFile->add_att("nnodes", "number of nodes" ) );
      assert( _ncFile->add_att("nfaces", "number of faces" ) );
      assert( _ncFile->add_att("ncells", "number of cells" ) );
 

      _nfaceRow  = _ncFile->add_dim("nface_row", nface_row);
      _nfaceCellsCol = _ncFile->add_dim("nfaceCells_col", nfaceCells_col);
      _nfaceNodesCol = _ncFile->add_dim("nfaceNodes_col", nfaceNodes_col);
      assert( _ncFile->add_att("nface_row", "row dimension of face connectivity" ) );
      assert( _ncFile->add_att("nfaceCells_col", "col dimension of faceCells connectivity" ) );
      assert( _ncFile->add_att("nfaceNodes_col", "col dimension of faceNodes connectivity" ) );

       if ( ninterface > 0 ){
           INTERFACE = true;
          _nInterface    = _ncFile->add_dim("nInterface", ninterface);
          assert( _ncFile->add_att("nInterface", "total interfaces") );
       }


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
     if ( BOUN_TYPE_DIM ){
        _boundarySize   = _ncFile->add_var("boundary_size", ncInt, _nBoun);
        _boundaryOffset = _ncFile->add_var("boundary_offset", ncInt, _nBoun);
        _boundaryID     = _ncFile->add_var("boundary_id", ncInt, _nBoun);
        _boundaryType   = _ncFile->add_var("boundary_type", ncChar, _nBoun, _charSize );
     }

    _interfaceGroup  = _ncFile->add_var("interface_group" , ncInt, _nmesh      );
    if ( NEIGH_MESH ){
       _interfaceSize   = _ncFile->add_var("interface_size"  , ncInt, _nNeighMesh );
       _interfaceOffset = _ncFile->add_var("interface_offset", ncInt, _nNeighMesh );
       _interfaceID     = _ncFile->add_var("interface_id"    , ncInt, _nNeighMesh ); 
    }

    _x  = _ncFile->add_var("x", ncDouble, _nnodes );
    _y  = _ncFile->add_var("y", ncDouble, _nnodes );
    _z  = _ncFile->add_var("z", ncDouble, _nnodes );

    _faceCellsRowCount = _ncFile->add_var("face_cells_row_count", ncInt, _nmesh );
    _faceCellsColCount = _ncFile->add_var("face_cells_col_count", ncInt, _nmesh );
    _faceNodesRowCount = _ncFile->add_var("face_nodes_row_count", ncInt, _nmesh );
    _faceNodesColCount = _ncFile->add_var("face_nodes_col_count", ncInt, _nmesh );

    _faceCellsRow = _ncFile->add_var("face_cells_row", ncInt, _nfaceRow );
    _faceCellsCol = _ncFile->add_var("face_cells_col", ncInt, _nfaceCellsCol );
    _faceNodesRow = _ncFile->add_var("face_nodes_row", ncInt, _nfaceRow );
    _faceNodesCol = _ncFile->add_var("face_nodes_col", ncInt, _nfaceNodesCol );
     
     if ( INTERFACE ){
        _fromIndices  = _ncFile->add_var("from_indices", ncInt, _nInterface);
        _toIndices    = _ncFile->add_var("to_indices"  , ncInt, _nInterface);
     }

    _bounBoolVar      = _ncFile->add_var("is_bounTypeDim_Valid", ncInt);
    _neighMeshBoolVar = _ncFile->add_var("is_neighMesh_Valid", ncInt);
    _interfaceBoolVar = _ncFile->add_var("is_interface_Valid", ncInt);



}

//assign values to NcVars
void 
NcDataWriter::set_var_values()
{
     //get variable values from meshes
     get_var_values();

     //adding attirbutes
     add_attributes();

     //write values
     write_values();


}



//getting values from meshes
//getting values from meshes
void 
NcDataWriter::get_var_values()
{

     assert( !_xVals );
     assert( !_yVals );
     assert( !_zVals );
     assert( !_fromIndicesVals );
     assert( !_toIndicesVals   );

    _xVals = new double [ _nnodes->size() ];
    _yVals = new double [ _nnodes->size() ];
    _zVals = new double [ _nnodes->size() ];

     if ( INTERFACE ){
       _fromIndicesVals  = new int [ _nInterface->size() ];
       _toIndicesVals    = new int [ _nInterface->size() ];
     }

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

        //boundary values
        get_boundary_vals( id );

       //interface values
        get_interface_vals( id );

       //x, y, z
        get_coords( id );

       //connectivities
       connectivities( id );

       //mappers
       if ( INTERFACE )
          mappers( id );

    }
}
//boundary face data
void
NcDataWriter::get_boundary_vals( int id )
{
       //boundary face
       const FaceGroupList& bounFaceList = _meshList.at(id)->getBoundaryFaceGroups();
       _boundaryGroupVals.push_back ( bounFaceList.size() );
       if ( BOUN_TYPE_DIM ) {
         for ( int boun = 0; boun < int(bounFaceList.size()); boun++ ) {
             _boundarySizeVals.push_back( bounFaceList.at(boun)->site.getCount() );
             _boundaryOffsetVals.push_back( bounFaceList.at(boun)->site.getOffset() );
             _boundaryIDVals.push_back( bounFaceList.at(boun)->id );
             _boundaryTypeVals.push_back( bounFaceList.at(boun)->groupType.c_str() );
              //assign values
             _boundaryType->set_cur(boun); 
              assert(  int(bounFaceList.at(boun)->groupType.size()) < MAX_CHAR );
             _boundaryType->put( _boundaryTypeVals[boun], 1, bounFaceList.at(boun)->groupType.size() );
          }
       }

}

//interface data
void
NcDataWriter::get_interface_vals( int id )
{
       //interface 
       const FaceGroupList& interfaceList = _meshList.at(id)->getInterfaceGroups();
       _interfaceGroupVals.push_back ( interfaceList.size() );
       if ( NEIGH_MESH ){
          for ( int interface = 0; interface < int( interfaceList.size() ); interface++ ){
             _interfaceSizeVals.push_back  (  interfaceList.at(interface)->site.getCount() );
             _interfaceOffsetVals.push_back(  interfaceList.at(interface)->site.getOffset() );
             _interfaceIDVals.push_back    (  interfaceList.at(interface)->id ); 
          }
       }

}

//coordinate values
void
NcDataWriter::get_coords( int id )
{
      int nn = _nodesCountVals.at(id);   
      const Mesh& mesh = *(_meshList.at(id));
      const Array<Mesh::VecD3>&  coord = mesh.getNodeCoordinates();
      for ( int n = 0; n < nn; n++ ){
         _xVals[n] = coord[n][0];
         _yVals[n] = coord[n][1];
         _zVals[n] = coord[n][2];
      }

       _x->put( _xVals, nn );
       _y->put( _yVals, nn );
       _z->put( _zVals, nn );

       _x->set_cur( nn );
       _y->set_cur( nn );
       _z->set_cur( nn );
}

//connectivities
void
NcDataWriter::connectivities( int id )
{
     //rows
     const Mesh& mesh = *(_meshList.at(id));
     const CRConnectivity& faceCells = mesh.getAllFaceCells();
     const CRConnectivity& faceNodes = mesh.getAllFaceNodes();

     int nRow = faceCells.getRow().getLength();
     _faceCellsRowCountVals.push_back( nRow );
     _faceNodesRowCountVals.push_back( nRow );
     _faceCellsRow->put( reinterpret_cast<int*> (faceCells.getRow().getData()), nRow );
     _faceNodesRow->put( reinterpret_cast<int*> (faceNodes.getRow().getData()), nRow );
     _faceCellsRow->set_cur( nRow );
     _faceNodesRow->set_cur( nRow );

    //cols
     int coldim = faceCells.getCol().getLength();
     _faceCellsColCountVals.push_back ( coldim );
     _faceCellsCol->put( reinterpret_cast<int*> (faceCells.getCol().getData()), coldim );
     _faceCellsCol->set_cur( coldim );

     //cols (faceNodes)
      coldim = faceNodes.getCol().getLength();
     _faceNodesColCountVals.push_back( coldim );
     _faceNodesCol->put( reinterpret_cast<int*> (faceNodes.getCol().getData()), coldim );
     _faceNodesCol->set_cur( coldim );

}


//MappersMap
void
NcDataWriter::mappers( int id )
{
     const StorageSite::MappersMap&   mappers = _meshList.at(id)->getCells().getMappers();
     const Mesh::GhostCellSiteMap&    ghostCellSiteMap = _meshList.at(id)->getGhostCellSiteMap();
     //loop over neighbour mesh to count interfaces
     StorageSite::MappersMap::const_iterator it_mapper;
     Mesh::GhostCellSiteMap::const_iterator  it;
 
     int  indx = mappers_index( id );
     
     for ( it = ghostCellSiteMap.begin(); it != ghostCellSiteMap.end(); it++ ){
         const StorageSite* site = it->second.get();
         it_mapper = mappers.find( site );
         int nend = it_mapper->second->getFromIndices().getLength();
         //int nend = it->second->getFromIndices().getLength();
   
         for ( int n = 0; n < nend; n++ ){
           _fromIndicesVals[indx] = it_mapper->second->getFromIndices()[n];  //or it->secon->toIndices.getLength() 
           _toIndicesVals[indx]   = it_mapper->second->getToIndices()[n];
            indx++;
         }
     }


}

//mappers index 
int
NcDataWriter::mappers_index( int mesh_end )
{
  int ninterface  = 0;
   for ( int id = 0; id < mesh_end; id++ ){
         const StorageSite::MappersMap&   mappers = _meshList.at(id)->getCells().getMappers();
         //loop over neighbour mesh to count interfaces
         StorageSite::MappersMap::const_iterator it;
         for ( it = mappers.begin(); it != mappers.end(); it++ )
             ninterface += it->second->getFromIndices().getLength();  //or it->secon->toIndices.getLength()

     }

    return ninterface;
}

//attributes
void 
NcDataWriter::add_attributes()
{
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
     if ( BOUN_TYPE_DIM ){
        assert( _boundarySize->add_att("boundary_size", " size of boundary" ) );
        assert( _boundaryOffset->add_att("boundary_offset", " offset of boundary" ) );
        assert( _boundaryID->add_att("boundary_id", " boundary id " ) );
        assert( _boundaryType->add_att("boundary_type", " type of boundary condition ") );
     }

     assert( _interfaceGroup->add_att("interface_group", " total interfaces") );

     if ( NEIGH_MESH ){
        assert( _interfaceSize->add_att("interface_size", " size of interface" ) );
        assert( _interfaceOffset->add_att("interface_offset", " offset of interface" ) );
        assert( _interfaceID->add_att("interface_id", " interface id " ) );
     }

     assert( _x->add_att("x", "x-coordinate") );
     assert( _y->add_att("y", "y-coordinate") );
     assert( _z->add_att("z", "z-coordinate") );

     assert( _faceCellsRowCount->add_att("face_cells_row_count", "count of row values of faceCells CRconnctivities") );
     assert( _faceNodesRowCount->add_att("face_nodes_row_count", "count of row values of faceNodes CRconnctivities") );
     assert( _faceCellsColCount->add_att("face_cells_col_count", "count of col values of faceCells CRconnctivities") );
     assert( _faceNodesColCount->add_att("face_nodes_col_count", "count of col values of faceNodes CRconnctivities") );

     assert( _faceCellsRow->add_att("face_cells_row", "row values of faceCells CRconnctivities") );
     assert( _faceNodesRow->add_att("face_nodes_row", "row values of faceNodes CRconnctivities") );
     assert( _faceCellsCol->add_att("face_cells_col", "col values of faceCells CRconnctivities") );
     assert( _faceNodesCol->add_att("face_nodes_col", "col values of faceNodes CRconnctivities") );
     if ( INTERFACE ){     
        assert( _fromIndices->add_att("from_indices", "trom indices from other neightbour mesh " ) );
        assert( _toIndices->add_att("to_indices",     "to  indices in current mesh") );
     }

}



//write values
void 
NcDataWriter::write_values()
{
    _dimension->put( &_dimensionVals[0], _nmesh->size()  );
    _meshID->put( &_meshIDVals[0], _nmesh->size() );
    _facesCount->put( &_facesCountVals[0], _nmesh->size() );
    _cellsCount->put( &_cellsCountVals[0], _nmesh->size() );
    _ghostCellsCount->put( &_ghostCellsCountVals[0], _nmesh->size() );
    _nodesCount->put( &_nodesCountVals[0], _nmesh->size() );
    _mapCount->put( &_mapCountVals[0], _nmesh->size() );
    _interiorFaceGroup->put( &_interiorFaceGroupVals[0], _nmesh->size() );

    _boundaryGroup->put(&_boundaryGroupVals[0], _nmesh->size() );
     if ( BOUN_TYPE_DIM ){
       _boundarySize->put( &_boundarySizeVals[0], _boundarySizeVals.size()  );
       _boundaryOffset->put( &_boundaryOffsetVals[0], _boundaryOffsetVals.size() );
       _boundaryID->put( &_boundaryIDVals[0], _boundaryIDVals.size() );
     }

    _interfaceGroup->put(&_interfaceGroupVals[0], _nmesh->size() );

     if ( NEIGH_MESH ){
       _interfaceSize->put( &_interfaceSizeVals[0], _interfaceSizeVals.size()  );
       _interfaceOffset->put( &_interfaceOffsetVals[0], _interfaceOffsetVals.size() );
       _interfaceID->put( &_interfaceIDVals[0], _interfaceIDVals.size() );
     }

     _faceCellsRowCount->put( &_faceCellsRowCountVals[0], _nmesh->size() );
     _faceCellsColCount->put( &_faceCellsColCountVals[0], _nmesh->size() );
     _faceNodesRowCount->put( &_faceNodesRowCountVals[0], _nmesh->size() );
     _faceNodesColCount->put( &_faceNodesColCountVals[0], _nmesh->size() );

    if ( INTERFACE ){
       _fromIndices->put( _fromIndicesVals, _nInterface->size() );
       _toIndices->put( _toIndicesVals, _nInterface->size() );
    }

     int boun_bool = int( BOUN_TYPE_DIM );
     int neigh_bool= int( NEIGH_MESH );
     int interface_bool = int ( INTERFACE );
    _bounBoolVar->put(  &boun_bool );
    _neighMeshBoolVar->put( &neigh_bool );
    _interfaceBoolVar->put( &interface_bool );


}
