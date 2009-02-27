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

#include "NcDataReader.h"
#include "netcdfcpp.h"
#include "mpi.h"


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



    if ( _ncFile        ) delete _ncFile;

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

  
}

//getting values
void
NcDataReader::get_var_values()
{

     allocate_vars();

    _dimension->get( _dimensionVals, _dimension->num_vals() );
    _meshID->get( _meshIDVals, _meshID->num_vals() ); 
    _facesCount->get( _facesCountVals, _facesCount->num_vals() );
    _cellsCount->get( _cellsCountVals, _cellsCount->num_vals() );
    _ghostCellsCount->get( _ghostCellsCountVals, _ghostCellsCount->num_vals() );
    _nodesCount->get( _nodesCountVals, _nodesCount->num_vals() );
    _mapCount->get ( _mapCountVals, _mapCount->num_vals() );
    _interiorFacesGroup->get( _interiorFacesGroupVals, _interiorFacesGroup->num_vals() );
      
     get_bndry_vals();

     //if ( MPI::COMM_WORLD.Get_rank() == 0) cout << _boundaryTypeVals[0] << "  " << _boundaryTypeVals[1] << endl;


}

void 
NcDataReader::allocate_vars()
{
    _dimensionVals = new  int[ _dimension->num_vals() ];
    _meshIDVals    = new  int[ _meshID->num_vals()    ];

    _facesCountVals = new int[ _facesCount->num_vals() ];
    _cellsCountVals = new int[ _cellsCount->num_vals() ];
    _ghostCellsCountVals = new int[ _ghostCellsCount->num_vals() ];
    _nodesCountVals  = new int[ _nodesCount->num_vals() ];
    _mapCountVals    = new int[ _mapCount->num_vals() ];
    _interiorFacesGroupVals = new int [ _interiorFacesGroup->num_vals() ];

    _boundaryGroupVals  = new int [ _boundaryGroup->num_vals()  ];
    _boundarySizeVals   = new int [ _boundarySize->num_vals()   ];
    _boundaryOffsetVals = new int [ _boundaryOffset->num_vals() ];
    _boundaryIDVals     = new int [ _boundaryID->num_vals()     ];
    _boundaryTypeVals.reserve(_nBoun);
    for ( int n = 0; n < _nBoun; n++)
         _boundaryTypeVals.push_back( new char[ _charSize ] );

}


void 
NcDataReader::get_bndry_vals()
{

   if ( _nBoun  > 0 ){
      _boundaryGroup->get (_boundaryGroupVals, _boundaryGroup->num_vals() );
      _boundarySize->get  (_boundarySizeVals,  _boundarySize->num_vals() );
      _boundaryOffset->get( _boundaryOffsetVals, _boundaryOffset->num_vals() );
      _boundaryID->get ( _boundaryIDVals, _boundaryID->num_vals() );
      for ( int n = 0; n < _nBoun; n++){
         _boundaryType->set_cur(n);
         _boundaryType->get( _boundaryTypeVals.at(n), 1, _charSize );
      }
   }
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


//         //then interface faces
//         for ( it_set = _interfaceSet.at(id).begin(); it_set != _interfaceSet.at(id).end(); it_set++){
//            int interfaceID = *it_set;
//            int size = int(_interfaceMap.at(id).count( interfaceID ) );
//            int offset = _interfaceOffsets.at(id)[ interfaceID];
//            _meshListLocal.at(id)->createInterfaceGroup( size, offset, interfaceID );
//            _meshListLocal.at(id)->createGhostCellSite( interfaceID, shared_ptr<StorageSite>( new StorageSite(size) ) );
//          }


        // mappers( id );


     }
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
       int boun_end = _boundaryGroupVals[id];
       cout << " proc_id = " << MPI::COMM_WORLD.Get_rank() <<    
       "  boun_end = " << boun_end << endl;
       for ( int boun = 0; boun < boun_end; boun++){
          int bndryID = _boundaryIDVals[boun];
          int size    = _boundarySizeVals[boun];
          int offset  = _boundaryOffsetVals[boun] ;
          string boundaryType (_boundaryTypeVals.at(boun) );
         _meshList.at(id)->createBoundaryFaceGroup( size, offset, bndryID, boundaryType);
       }
 
}

// void 
// NcDataReader::mappers( int id )
// {
// 
// }









