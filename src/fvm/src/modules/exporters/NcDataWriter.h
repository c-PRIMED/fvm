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

private :
    NcDataWriter( const NcDataWriter& nc_writer);

    void  init();
    void  setNcFile();
    void  setDims();
    void  setVars();
    void  set_var_values();

     const MeshList& _meshList;
     string  _fname;

     NcFile   *_ncFile;
     //NcDims
     NcDim    *_nmesh;
     NcDim    *_nBoun;
     NcDim    *_charSize;
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


     const int MAX_CHAR;

};

#endif
