//
// C++ Interface: NcDataREADER
//
// Description: 
//
//
// Author: yildirim,,, <yildirim@prism>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef NCDATAREADER_H
#define NCDATAREADER_H

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


class NcDataReader {

public :

    NcDataReader( const string& fname );
    
    void    read();

    ~NcDataReader();

private :

    void  init();
    void  setNcFile();
    void  getDims();
    void  getVars();
    void  get_var_values();

    void  meshList();
    void  mappers( int id );

    void storage_sites( int id );
    void boundary_faces( int id );
    void allocate_vars();
    void get_bndry_vals();



    string _fname;


    //netcdf variables    
    NcFile   *_ncFile;
    int  _nmesh;
    int  _nBoun;
    int  _charSize;


    NcVar *_dimension;
    NcVar *_meshID;
    NcVar *_facesCount;
    NcVar *_cellsCount;
    NcVar *_ghostCellsCount;
    NcVar *_nodesCount;
    NcVar *_mapCount;
    NcVar *_interiorFacesGroup;

    NcVar*  _boundaryGroup;
    NcVar*  _boundarySize;
    NcVar*  _boundaryOffset;
    NcVar*  _boundaryID;
    NcVar*  _boundaryType;

    int*  _dimensionVals;
    int*  _meshIDVals;

    int *_facesCountVals;
    int *_cellsCountVals;
    int *_ghostCellsCountVals;
    int *_nodesCountVals;
    int *_mapCountVals;
    int *_interiorFacesGroupVals;


    int *_boundaryGroupVals;
    int *_boundarySizeVals;
    int *_boundaryOffsetVals;
    int * _boundaryIDVals;
    vector<char*>  _boundaryTypeVals;


    MeshList   _meshList;
};

#endif
