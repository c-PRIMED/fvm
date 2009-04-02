//
// C++ Interface: NcDataWriter
//
// Description: 
//
//
// Author: yildirim,,, <yildirim@prism>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef NCDATAWRITER_H
#define NCDATAWRITER_H

#include <string>
#include <vector>
#include "Mesh.h"
using namespace std;

/**
	@author yildirim,,, <yildirim@prism>
*/

class NcFile;
class NcDim;
class NcVar;

class NcDataWriter {

public :

    NcDataWriter(const MeshList& meshes, const string& fname);
    ~NcDataWriter();

    void  record();

     static int _writeAction;

private :
    NcDataWriter( const NcDataWriter& nc_writer);

    void  init();
    void  setNcFile();
    void  setDims();
    void  setVars();
    void  set_var_values();

    void  get_var_values();
    void  get_boundary_vals( int id );
    void  get_interface_vals( int id );
    void  get_coords( int id );
    void  connectivities( int id );
    void  mappers( int id );
    int   mappers_index( int mesh_end );

    void  add_attributes();
    void  write_values();

    




     const MeshList _meshList;
     string  _fname;

     NcFile   *_ncFile;
     //NcDims
     NcDim    *_nmesh;
     NcDim    *_nBoun;
     NcDim    *_charSize;
     NcDim    *_nNeighMesh;
     NcDim    *_nnodes;
     NcDim    *_nfaces;
     NcDim    *_ncells;

     NcDim    *_nfaceRow;
     NcDim    *_nfaceCellsCol;
     NcDim    *_nfaceNodesCol;
    
     NcDim    *_nInterface;

     //NcVars
     NcVar*  _dimension;
     NcVar*  _meshID;
     NcVar*  _facesCount;
     NcVar*  _cellsCount;
     NcVar*  _ghostCellsCount;
     NcVar*  _nodesCount;
     NcVar*  _mapCount;
     NcVar*  _interiorFaceGroup;
     
     NcVar*  _boundaryGroup;
     NcVar*  _boundarySize;
     NcVar*  _boundaryOffset;
     NcVar*  _boundaryID;
     NcVar*  _boundaryType;

     NcVar*  _interfaceGroup;
     NcVar*  _interfaceSize;
     NcVar*  _interfaceOffset;
     NcVar*  _interfaceID;

     NcVar*  _x;
     NcVar*  _y;
     NcVar*  _z;
    
     NcVar* _faceCellsRowCount;
     NcVar* _faceCellsColCount;
     NcVar* _faceNodesRowCount;
     NcVar* _faceNodesColCount;

     NcVar* _faceCellsRow;
     NcVar* _faceCellsCol;
     NcVar* _faceNodesRow;
     NcVar* _faceNodesCol;

     NcVar* _gatherIndices;
     NcVar* _scatterIndices;


     NcVar* _bounBoolVar;
     NcVar* _neighMeshBoolVar;
     NcVar* _interfaceBoolVar;

     //variable values
     vector< int >  _dimensionVals;
     vector< int >  _meshIDVals;

     vector< int >  _facesCountVals;
     vector< int >  _cellsCountVals;
     vector< int >  _ghostCellsCountVals;
     vector< int >  _nodesCountVals;
     vector< int >  _mapCountVals;
     vector< int > _interiorFaceGroupVals;

     vector< int > _boundaryGroupVals;
     vector< int > _boundarySizeVals;
     vector< int > _boundaryOffsetVals;
     vector< int > _boundaryIDVals;
     vector< const char* > _boundaryTypeVals;


     vector< int > _interfaceGroupVals;
     vector< int > _interfaceSizeVals;
     vector< int > _interfaceOffsetVals;
     vector< int > _interfaceIDVals;

     double  *_xVals; 
     double  *_yVals;
     double  *_zVals;

     vector< int > _faceCellsRowCountVals;
     vector< int > _faceCellsColCountVals;
     vector< int > _faceNodesRowCountVals;
     vector< int > _faceNodesColCountVals;

     int* _gatherIndicesVals;
     int* _scatterIndicesVals;

     const int MAX_CHAR;
     bool BOUN_TYPE_DIM;
     bool NEIGH_MESH;
     bool INTERFACE;

};

#endif
