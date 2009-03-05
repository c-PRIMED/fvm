//
// C++ Implementation: NcDataReader
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
#include <algorithm>
#include <numeric>

#include "NcDataReader.h"
#include "netcdfcpp.h"
#include "mpi.h"
#include "OneToOneIndexMap.h"


NcDataReader::NcDataReader( const string& fname )
:_fname( fname )
{
    init();
}


NcDataReader::~NcDataReader()
{

    if ( _dimensionVals ) delete [] _dimensionVals;
    if ( _meshIDVals    ) delete [] _meshIDVals;

    if ( _facesCountVals ) delete []  _facesCountVals;
    if ( _cellsCountVals ) delete [] _cellsCountVals;
    if ( _ghostCellsCountVals ) delete [] _ghostCellsCountVals;
    if (_nodesCountVals ) delete [] _nodesCountVals;
    if ( _mapCountVals  ) delete [] _mapCountVals;
    if ( _interiorFacesGroupVals ) delete [] _interiorFacesGroupVals;

    if ( _boundaryGroupVals )  delete [] _boundaryGroupVals;
    if ( _boundarySizeVals )   delete [] _boundarySizeVals;
    if ( _boundaryOffsetVals ) delete [] _boundaryOffsetVals;
    if ( _boundaryIDVals )     delete [] _boundaryIDVals;

    vector< char* >::iterator it_char;
    for ( it_char = _boundaryTypeVals.begin(); it_char != _boundaryTypeVals.end(); it_char++ )
         delete [] *it_char;


    if ( _interfaceGroupVals )  delete [] _interfaceGroupVals;
    if ( _interfaceSizeVals )   delete [] _interfaceSizeVals;
    if ( _interfaceOffsetVals ) delete [] _interfaceOffsetVals;
    if ( _interfaceIDVals )     delete [] _interfaceIDVals;
  
    if ( _xVals )   delete [] _xVals;
    if ( _yVals )   delete [] _yVals;
    if ( _zVals )   delete [] _zVals;

    if ( _faceCellsRowVals ) delete [] _faceCellsRowVals;
    if ( _faceNodesRowVals ) delete [] _faceNodesRowVals;
    if ( _faceCellsColVals ) delete [] _faceCellsColVals;
    if ( _faceNodesColVals ) delete [] _faceNodesColVals;


    if ( _fromIndicesVals  ) delete [] _fromIndicesVals;
    if ( _toIndicesVals    ) delete [] _toIndicesVals;
   

    if ( _ncFile        ) delete _ncFile;   

    MeshList::iterator it_mesh;
   for ( it_mesh = _meshList.begin(); it_mesh != _meshList.end(); it_mesh++)
       delete  *it_mesh;


}



void
NcDataReader::read()
{
   setNcFile();
   getDims();
   getVars();
   get_var_values();

   meshList();

}


                     //PRIVATE FUNCTIONS
void
NcDataReader::init()
{
    _ncFile = NULL;
    _dimensionVals = NULL;
    _meshIDVals  = NULL;
    _facesCountVals = NULL;
    _cellsCountVals = NULL;
    _ghostCellsCountVals = NULL;
    _nodesCountVals = NULL;
    _mapCountVals = NULL;
    _interiorFacesGroupVals = NULL;

    _boundaryGroupVals = NULL;
    _boundarySizeVals  = NULL;
    _boundaryOffsetVals= NULL;
    _boundaryIDVals    = NULL;

    _interfaceGroupVals  = NULL;
    _interfaceSizeVals   = NULL;
    _interfaceOffsetVals = NULL;
    _interfaceIDVals     = NULL;

    _xVals = NULL;
    _yVals = NULL;
    _zVals = NULL;

    _faceCellsRowCountVals = NULL;
    _faceNodesRowCountVals = NULL;
    _faceCellsColCountVals = NULL;
    _faceNodesColCountVals = NULL;

    _faceCellsRowVals = NULL;
    _faceNodesRowVals = NULL;
    _faceCellsColVals = NULL;
    _faceNodesColVals = NULL;

    _fromIndicesVals = NULL;
    _toIndicesVals   = NULL;
}		   

//Setting NcFile
void
NcDataReader::setNcFile()
{
   assert( !_ncFile );
   _ncFile = new NcFile( _fname.c_str(), NcFile::ReadOnly );
   assert ( _ncFile->is_valid() );

}  


//Getting dimension from NcFile
void
NcDataReader::getDims()
{
   _nmesh    = _ncFile->get_dim("nmesh")->size();
   _nBoun    = _ncFile->get_dim("boun_type_dim")->size();
   _charSize = _ncFile->get_dim("char_dim")->size();
   _nNeighMesh = _ncFile->get_dim("nNeighMesh")->size();
   _nnodes     = _ncFile->get_dim("nnodes")->size();
   _nfaceRow   = _ncFile->get_dim("nface_row")->size();
   _nfaceCellsCol = _ncFile->get_dim("nfaceCells_col")->size();
   _nfaceNodesCol = _ncFile->get_dim("nfaceNodes_col")->size();
   _nInterface    = _ncFile->get_dim("nInterface")->size();


}

//getting NCVar-s
void 
NcDataReader::getVars()
{
    _dimension =  _ncFile->get_var("dimension");
    _meshID    =  _ncFile->get_var("mesh_id"  ); 
    _facesCount = _ncFile->get_var("faces_count");
    _cellsCount = _ncFile->get_var("cells_count");
    _ghostCellsCount = _ncFile->get_var("ghost_cells_count");
    _nodesCount = _ncFile->get_var("nodes_count");
    _mapCount   = _ncFile->get_var("map_count");
    _interiorFacesGroup = _ncFile->get_var("interior_faces_group");

    _boundaryGroup   = _ncFile->get_var("boundary_group");
    _boundarySize    = _ncFile->get_var("boundary_size");
    _boundaryOffset  = _ncFile->get_var("boundary_offset");
    _boundaryID      = _ncFile->get_var("boundary_id"); 
    _boundaryType    = _ncFile->get_var("boundary_type"); 


    _interfaceGroup   = _ncFile->get_var("interface_group");
    _interfaceSize    = _ncFile->get_var("interface_size");
    _interfaceOffset  = _ncFile->get_var("interface_offset");
    _interfaceID      = _ncFile->get_var("interface_id"); 

    _x = _ncFile->get_var("x");
    _y = _ncFile->get_var("y");
    _z = _ncFile->get_var("z");     

    _faceCellsRowCount = _ncFile->get_var("face_cells_row_count");
    _faceCellsColCount = _ncFile->get_var("face_cells_col_count");
    _faceNodesRowCount = _ncFile->get_var("face_nodes_row_count");
    _faceNodesColCount = _ncFile->get_var("face_nodes_col_count");


    _faceCellsRow = _ncFile->get_var("face_cells_row");
    _faceCellsCol = _ncFile->get_var("face_cells_col");
    _faceNodesRow = _ncFile->get_var("face_nodes_row");
    _faceNodesCol = _ncFile->get_var("face_nodes_col");

    _fromIndices = _ncFile->get_var("from_indices");
    _toIndices   = _ncFile->get_var("to_indices");


  
}

//getting values
void
NcDataReader::get_var_values()
{

     allocate_vars();

    _dimension->get( _dimensionVals, _nmesh );
    _meshID->get( _meshIDVals      , _nmesh ); 
    _facesCount->get( _facesCountVals, _nmesh );
    _cellsCount->get( _cellsCountVals, _nmesh );
    _ghostCellsCount->get( _ghostCellsCountVals, _nmesh );
    _nodesCount->get( _nodesCountVals, _nmesh );
    _mapCount->get ( _mapCountVals   , _nmesh );
    _interiorFacesGroup->get( _interiorFacesGroupVals, _nmesh );
      
     get_bndry_vals();
     get_interface_vals(); 
     get_coord_vals();
     get_connectivity_vals();
     get_mapper_vals();

}

void 
NcDataReader::allocate_vars()
{
    _dimensionVals = new  int[ _nmesh ];
    _meshIDVals    = new  int[ _nmesh    ];

    _facesCountVals = new int[ _nmesh ];
    _cellsCountVals = new int[ _nmesh ];
    _ghostCellsCountVals = new int[ _nmesh ];
    _nodesCountVals  = new int[ _nmesh ];
    _mapCountVals    = new int[ _nmesh ];
    _interiorFacesGroupVals = new int [ _nmesh ];

    _boundaryGroupVals  = new int [ _nmesh ];
    _boundarySizeVals   = new int [ _nBoun ];
    _boundaryOffsetVals = new int [ _nBoun ];
    _boundaryIDVals     = new int [ _nBoun ];
    _boundaryTypeVals.reserve(_nBoun);
    for ( int n = 0; n < _nBoun; n++)
         _boundaryTypeVals.push_back( new char[ _charSize ] );

    _interfaceGroupVals  = new int [ _nmesh  ];
    _interfaceSizeVals   = new int [ _nNeighMesh ];
    _interfaceOffsetVals = new int [ _nNeighMesh ];
    _interfaceIDVals     = new int [ _nNeighMesh ];

    _xVals = new double [ _nnodes ];
    _yVals = new double [ _nnodes ];
    _zVals = new double [ _nnodes ];
 


   _faceCellsRowCountVals = new int [ _nmesh ];
   _faceCellsColCountVals = new int [ _nmesh ];
   _faceNodesRowCountVals = new int [ _nmesh ];
   _faceNodesColCountVals = new int [ _nmesh ];

   _faceCellsRowVals = new int [ _nfaceRow ];
   _faceCellsColVals = new int [ _nfaceCellsCol ];
   _faceNodesRowVals = new int [ _nfaceRow ];
   _faceNodesColVals = new int [ _nfaceNodesCol ];
   
   _fromIndicesVals  = new int [ _nInterface ];
   _toIndicesVals    = new int [ _nInterface ];

}

//get boundary values
void 
NcDataReader::get_bndry_vals()
{

   if ( _nBoun  > 0 ){
      _boundaryGroup->get (_boundaryGroupVals, _nmesh );
      _boundarySize->get  (_boundarySizeVals,  _nBoun );
      _boundaryOffset->get( _boundaryOffsetVals, _nBoun );
      _boundaryID->get ( _boundaryIDVals, _nBoun );
      for ( int n = 0; n < _nBoun; n++){
         _boundaryType->set_cur(n);
         _boundaryType->get( _boundaryTypeVals.at(n), 1, _charSize );
      }
   }
}



void 
NcDataReader::get_interface_vals()
{
   if ( _nNeighMesh  > 0 ){
      _interfaceGroup->get (_interfaceGroupVals, _nmesh );
      _interfaceSize->get  (_interfaceSizeVals,  _nNeighMesh );
      _interfaceOffset->get( _interfaceOffsetVals, _nNeighMesh );
      _interfaceID->get ( _interfaceIDVals, _nNeighMesh );
   }


}

//get coordinates of nodes
void
NcDataReader::get_coord_vals()
{

    _x->get( _xVals, _nnodes );
    _y->get( _yVals, _nnodes );
    _z->get( _zVals, _nnodes );


}

//get coordinates of nodes
void
NcDataReader::get_connectivity_vals()
{
    _faceCellsRowCount->get( _faceCellsRowCountVals, _nmesh );
    _faceNodesRowCount->get( _faceNodesRowCountVals, _nmesh );
    _faceCellsColCount->get( _faceCellsColCountVals, _nmesh );
    _faceNodesColCount->get( _faceNodesColCountVals, _nmesh );

    _faceCellsRow->get( _faceCellsRowVals, _nfaceRow );
    _faceNodesRow->get( _faceNodesRowVals, _nfaceRow );
    _faceCellsCol->get( _faceCellsColVals, _nfaceCellsCol );
    _faceNodesCol->get( _faceNodesColVals, _nfaceNodesCol );


}

//get mapper values
void
NcDataReader::get_mapper_vals()
{
    _fromIndices->get( _fromIndicesVals, _nInterface );
      _toIndices->get( _toIndicesVals  , _nInterface );

}



//forming MeshList
void
NcDataReader::meshList()
{

     for ( int id = 0; id < _nmesh; id++ ){

        _meshList.push_back(  new Mesh( _dimensionVals[id], _meshIDVals[id] )  );
         //storage sites
         storage_sites( id);
        //interior faces
         _meshList.at(id)->createInteriorFaceGroup( _interiorFacesGroupVals[id] );
        //boundary faces
        if ( _nBoun > 0 )
           boundary_faces( id );

        if ( _nNeighMesh > 0 )
           interfaces( id );

        coords( id );
        face_cells( id );
        face_nodes( id );

     }

    mappers();
}


void 
NcDataReader::storage_sites( int id )
{
         StorageSite& faces = _meshList.at(id)->getFaces();
         StorageSite& cells = _meshList.at(id)->getCells();
         StorageSite& nodes = _meshList.at(id)->getNodes();

         faces.setCount( _facesCountVals[id] );
         cells.setCount( _cellsCountVals[id], _ghostCellsCountVals[id] );
         nodes.setCount( _nodesCountVals[id] );
}



void
NcDataReader::boundary_faces( int id )
{
       //boundary faces
       int indx = accumulate( _boundaryGroupVals, _boundaryGroupVals+id, 0 );
       int nboun   = _boundaryGroupVals[id];
       for ( int boun = 0; boun < nboun; boun++){
          int bndryID = _boundaryIDVals[indx];
          int size    = _boundarySizeVals[indx];
          int offset  = _boundaryOffsetVals[indx] ;
          string boundaryType (_boundaryTypeVals.at(indx) );
         _meshList.at(id)->createBoundaryFaceGroup( size, offset, bndryID, boundaryType);
         indx++;
       }
 
}

void
NcDataReader::interfaces( int id )
{
     //then interface faces
     int indx = accumulate( _interfaceGroupVals, _interfaceGroupVals+id,0);
     int  ninterfaces = _interfaceGroupVals[id];
     for ( int interface = 0; interface < ninterfaces; interface++ ){
         int interfaceID = _interfaceIDVals[indx];
         int size        = _interfaceSizeVals[indx]; 
         int offset      = _interfaceOffsetVals[indx];
         _meshList.at(id)->createInterfaceGroup( size, offset, interfaceID );
         _meshList.at(id)->createGhostCellSite( interfaceID, shared_ptr<StorageSite>( new StorageSite(size) ) );
         indx++;
     }



}

void
NcDataReader::coords( int id )
{
     int nnodes = _nodesCountVals[id]; 
     int indx =  accumulate(_nodesCountVals, _nodesCountVals+id,0);
     shared_ptr< Array<Mesh::VecD3> >  coord( new Array<Mesh::VecD3>( nnodes ) );

     for ( int n = 0; n < nnodes; n++ ){
          (*coord)[n][0] = _xVals[indx];
          (*coord)[n][1] = _yVals[indx];
          (*coord)[n][2] = _zVals[indx];
          indx++;
      }

       _meshList.at(id)->setCoordinates( coord );


}

//connectivities
void
NcDataReader::face_cells( int id )
{
     //faceCells
    int nfaces = _facesCountVals[id];
     StorageSitePtr  rowSite( new StorageSite( nfaces ) );
     StorageSitePtr  colSite( new StorageSite( _cellsCountVals[id], _ghostCellsCountVals[id] ) );

     CRConnectivityPtr  faceCells ( new CRConnectivity( *rowSite, *colSite ) );

     faceCells->initCount();

     for ( int n = 0; n < nfaces; n++ ) 
         faceCells->addCount(n,2);  // two cells around a face

     faceCells->finishCount();

     int indx =  accumulate(_faceCellsColCountVals, _faceCellsColCountVals+id, 0);
     for ( int n = 0; n < nfaces; n++ ){
         for ( int cell = 0; cell < 2; cell++){ //two cells around a face always
             faceCells->add( n, _faceCellsColVals[indx] );
         indx++;
         }
     }

     faceCells->finishAdd(); 
    _meshList.at(id)->setFaceCells( faceCells );
}

//connectivities
void
NcDataReader::face_nodes( int id )
{
     //faceNodes
    int nfaces = _facesCountVals[id];
     StorageSitePtr  rowSite( new StorageSite( nfaces ) );
     StorageSitePtr  colSite( new StorageSite( _nodesCountVals[id]) );

     CRConnectivityPtr  faceNodes ( new CRConnectivity( *rowSite, *colSite ) );

     faceNodes->initCount();
     int node_count = _faceNodesRowVals[1] - _faceNodesRowVals[0];   
     for ( int n = 0; n < nfaces; n++ ) 
         faceNodes->addCount(n, node_count);

     faceNodes->finishCount();

     int indx =  accumulate(_faceNodesColCountVals, _faceNodesColCountVals+id, 0);
     for ( int n = 0; n < nfaces; n++ ){
         for ( int node = 0; node < node_count; node++){
            faceNodes->add( n, _faceNodesColVals[indx] );
            indx++;
         }
     }

     faceNodes->finishAdd(); 
     _meshList.at(id)->setFaceNodes  ( faceNodes );


}



void 
NcDataReader::mappers(  )
{

    int indx = 0;
    for ( int id = 0; id < _nmesh; id++){
        StorageSite::MappersMap&   cellMappers = _meshList.at(id)->getCells().getMappers();
       //loop over mesh interfaces
        int offset = accumulate( _interfaceGroupVals, _interfaceGroupVals+id,0);
        for ( int n = 0; n  < _interfaceGroupVals[id]; n++){
           int neighMeshID =  _interfaceIDVals[ offset + n ];
           int size = _interfaceSizeVals[ offset + n ];
           ArrayIntPtr  fromIndices( new Array<int>( size ) );
           ArrayIntPtr  toIndices  ( new Array<int>( size ) );


          //get portion values
          for ( int i = 0; i < size; i++){
              (*fromIndices)[i] = _fromIndicesVals[indx];
              (*toIndices)[i]   = _toIndicesVals[indx];
              indx++;
           }
           shared_ptr< OneToOneIndexMap >  oneToOneMapPtr( new OneToOneIndexMap(fromIndices, toIndices)  );
           cellMappers[ _meshList.at(id)->getGhostCellSite( neighMeshID ) ] =  oneToOneMapPtr;
       }
    }


}










