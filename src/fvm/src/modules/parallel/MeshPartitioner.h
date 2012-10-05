// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

//
// C++ Interface: MeshPartitioner
//
// Description: 
//
//
// Author: yildirim,,, <yildirim@prism>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MeshPartitioner_H
#define MeshPartitioner_H

#include <vector>
#include <map>
#include <string>
#include "Mesh.h"
#include "Array.h"
#include "Vector.h"
#include <mpi.h>
#include<set>
#include<fstream>

using namespace std;

/**
	@author yildirim,,, <yildirim@cfm>
*/

//Warning globalToLocal seems have a bug..........


class MeshPartitioner{

public:
    //typedef map<const StorageSite*, shared_ptr<OneToOneIndexMap> > MappersMap;
    typedef   shared_ptr< StorageSite >  StorageSitePtr; 
    typedef   shared_ptr< CRConnectivity > CRConnectivityPtr;
    typedef   shared_ptr< Array<int> >     ArrayIntPtr;
    typedef   shared_ptr< Array<Mesh::VecD3> >   ArrayVecD3Ptr;
    typedef   shared_ptr< Mesh >  MeshPtr;
    typedef   map<int, int>       IntMap;
    typedef   map<int, vector<int> > VectorMap;
    typedef   pair<int,int>       IntPair;

    enum ETYPE{ TRI = 1, QUAD = 2, TETRA = 3, HEXA = 4 };
    enum WTYPE{ NOWEIGHTS = 0, WEIGHTS_ONLY_EDGES = 1, WEIGTHS_ONLY_VERTICES  = 2,
                WEIGHTS_BOTH_VERTICES_EDGES = 3};
    enum NUMFLAG{ C_STYLE = 0, FORTRAN_STYLE = 1 };
    enum CELLTYPE{ INTERIOR = 1, GHOST_BOUNDARY_CELL = 2, GHOST_INTERFACE_CELL};
    enum PARTTYPE{ PARMETIS = 1, FIEDLER =2};

   explicit MeshPartitioner( const MeshList& mesh_list, vector<int> npart, vector<int> eType);
   ~MeshPartitioner();

    void partition();
    void fiedler_order( const string& fname );

    void mesh();
    const MeshList&  meshList() const { return _meshListLocal; }
    void  extractBoundaryMesh();
    const Mesh& getBoundaryMesh() const { return *_bMesh; }


    void debug_print();
#if 0
    void dumpTecplot();
#endif    
    void mesh_xdmfplot();

    // set property methods
    void setWeightType( int weight_type );
    void setNumFlag( int num_flag);

    //clean up memory 
    void cleanup();
    void isCleanup(bool clean_up) { _cleanup = clean_up;  }
    void isDebug  (bool debug   ) { _cleanup = false;
                                    _debugMode = debug;     }
    //this need to be public since it will be called by python script
    void DEBUG_compute_elem_dist();
    void DEBUG_parmetis_mesh();
    void DEBUG_fiedler_partition();
    void DEBUG_elem_connectivity();
    void DEBUG_map_part_elms();
    void DEBUG_count_elems_part();
    void DEBUG_exchange_part_elems();
    void DEBUG_mapBounIDAndCell();
    void DEBUG_resize_elem();
    void DEBUG_CRConnectivity_cellParts();
    void DEBUG_CRConnectivity_faceParts();
    void DEBUG_interfaces();
    void DEBUG_faceCells_faceNodes();
    void DEBUG_non_interior_cells();
    void DEBUG_preserve_cell_order();
    void DEBUG_order_faceCells_faceNodes();
    void DEBUG_coordinates();
    void DEBUG_exchange_interface_meshes();
    void DEBUG_mesh();
    void DEBUG_local_global();
    void DEBUG_cellcells_global();
    void DEBUG_globalCellID_procID_map();
    void DEBUG_gatherCellsLevel1_partID_map();
    void DEBUG_level1_scatter_gather_cells();
    void DEBUG_CRConnectivity_cellCells2();

private:

   MeshPartitioner( const MeshPartitioner& part_mesh ); //dont allow copy constructor
    
   void init();
   void compute_elem_dist();
   void elem_connectivity();
   void parmetis_mesh();
   void fiedler_partition();
   int  get_local_nodes(int id);
   void set_eptr_eind(int id);
   void map_part_elms();
   void count_elems_part();
   void exchange_part_elems();
   void shift_sum_row();
   void mesh_setup();
   int global_offset();
   void set_local_global();
   void set_cellcells_global();
   void cellcells_global_extension();
   void globalCellID_procID_map();
   void gatherCellsLevel1_partID_map();
   void level1_scatter_gather_cells();


   void CRConnectivity_cellParts();
   void mapBounIDAndCell(int id);
   void resize_elem(int id);
   void CRConnectivity_faceParts(); 
   void faceCells_faceNodes(); 

   void interfaces();
   void coordinates();
   int  count_interior_faces( int id );
   void order_faceCells_faceNodes(); 
   void non_interior_cells();
   void exchange_interface_meshes();
   void mappers();
   void create_window( int id );
   void free_window( );
   void fence_window();
   int  get_window_displ( int id, int neigh_mesh_id );
   void construct_mesh( int id );
   void setMeshColors();
   void preserve_cell_order();
   void scatterCellsLevel1();



   void mesh_xdmf_header();
   void mesh_file();
   void mesh_tecplot();


   void cleanup_follow_exchange_part_elems();
   void cleanup_follow_faceCells_faceNodes();
   void cleanup_follow_mappers();

   void CRConnectivityPrint( const CRConnectivity& conn, int procID, const string& name );
   void CRConnectivityPrintFile( const CRConnectivity& conn, const string& name );



   void  debug_file_open( const string& fname );
   void  debug_file_close();

   const MeshList _meshList;
   vector<int> _nPart;
   vector<int> _totElems;
   vector<int> _totElemsAndGhosts;
   vector< Array<int>* > _elemDist;   //store dist element number for each Mesh
   vector< Array<int>* > _globalIndx; //store where partition starts in glbal array
   vector< int* > _ePtr;  
   vector< int* > _eInd;
   vector< int* > _eElm;
   vector< int* > _elmWght;
   vector< int  > _wghtFlag;
   vector< int  > _numFlag;
   vector< int  > _ncon;
   vector< int  > _ncommonNodes;
   vector< int > _eType;
   vector< float* > _tpwgts;
   vector< float* > _ubvec;
   int _options;
   int _procID;
   //output variables
   vector< int > _edgecut;
   vector< int* > _part;

  //
   vector< multimap<int,int> > _mapPartAndElms;
   vector< int* > _row;
   vector< int* > _col;
   vector< int* > _elem;           //just interior (local)
   vector< set<int> > _elemSet;   
   vector< int* > _elemWithGhosts; //interior+boundary (local)
   

   vector<int> _nelems;
   vector<int> _nelemsWithGhosts;
   vector<int> _colDim;
   int _nmesh;

   vector< StorageSitePtr > _cellSiteGlobal;
   vector< StorageSitePtr > _cellSite;
   vector< StorageSitePtr > _nodeSite;
   vector< StorageSitePtr > _faceSite;
   vector< StorageSitePtr > _partSite;


   vector< const CRConnectivity* >  _faceCellsGlobal;      //global mesh
   vector< const CRConnectivity* >  _faceNodesGlobal;      //global mesh

   vector< CRConnectivityPtr >  _faceCells; //belongs to this partition
   vector< CRConnectivityPtr >  _faceNodes; //belongs to this partition
   vector< CRConnectivityPtr >  _cellCells; //belongs to this partition
   vector< CRConnectivityPtr >  _cellNodes; //belongs to this partition
   vector< CRConnectivityPtr >  _faceCellsOrdered;  
   vector< CRConnectivityPtr >  _faceNodesOrdered;

   vector< CRConnectivityPtr >  _cellParts;
   vector< CRConnectivityPtr >  _partCells;
   vector< CRConnectivityPtr >  _faceParts; 
   vector< CRConnectivityPtr >  _partFaces; //transpose of _faceParts
   vector< CRConnectivityPtr >  _partNodes;


   vector< set<int>  >     _boundarySet;
   vector< multimap<int,int> > _mapBounIDAndCell;  //belongs to this partition (global numbering)
   vector< map<int,string> >   _mapBounIDAndBounType; //belongs to this partition

   vector< ArrayVecD3Ptr > _coord;     //belongs to this partition
   vector< set<int> > _interfaceSet;
   vector< multimap<int,int>  > _interfaceMap;   //belongs to this partition (local numbering)

   vector< set<int> > _nonInteriorCells;      //local numbering

   vector< map<int,int> > _bndryOffsets;
   vector< map<int,int> > _interfaceOffsets;
   vector< vector<int> > _cellToOrderedCell;
   vector< multimap<int,int> > _globalToLocalMappers; 
   vector< map<int,int> > _localToGlobalMappers;
   map<int,int>           _globalToLocal; 
   map<int,int>           _cellToPreservedOrderCell;
   map<int,int>           _cellsLevel1PartID;
   map<int,int>           _gatherCellsLevel1PartIDMap;

  //variables aglomorated for MPI communications
   vector< ArrayIntPtr >  _interfaceMeshCounts; //ArrayIntPtr[0] = value , 0th mesh has total "value" neigh interfaced meshes
   vector< ArrayIntPtr >  _offsetInterfaceCells; //this start location of interfaces, need at communuctation a lot 
   vector< ArrayIntPtr >  _procTotalInterfaces; //store total interfaces for each process
   vector< ArrayIntPtr >  _interfaceMeshIDs;     //
   vector< ArrayIntPtr >  _ghostCellsGlobal;   //stored in contigous memory (global numbering) local data
   vector< ArrayIntPtr >  _ghostCellsLocal;      //local numbering local data
   vector< int > _windowSize;


   vector< vector< ArrayIntPtr > > _fromIndices;
   vector< vector< ArrayIntPtr > > _toIndices;
   ArrayIntPtr  _fiedlerMap; // map original cell number to  fiedler numbering
   PARTTYPE _partTYPE;

   MPI::Win  _winGlobal;
   MPI::Win  _winLocal;
   bool  _cleanup;
   bool  _debugMode;
   ofstream   _debugFile;
   vector< Mesh* >  _meshListLocal;  //all is for you
   Mesh*  _bMesh; //this is boundary mesh of whole object
};

#endif
