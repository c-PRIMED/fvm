// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CELLMARK_H_
#define _CELLMARK_H_


#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "StorageSite.h"
#include "FlowFields.h"
#include "GeomFields.h"
#include "Vector.h"
#include "NumType.h"
#include "Octree.h"
#include "CRConnectivity.h"
#include "Array.h"

typedef Vector<double, 3> VecD3;
typedef Array<VecD3> VecD3Array;

int inCell(const int cellIndex, 
	   const Vector<double, 3>& point, 
	   const CRConnectivity& faceCells,
	   const CRConnectivity& cellFaces,
	   const VecD3Array& faceArea,
	   const VecD3Array& faceCentroid);


void reportCellMark (const Mesh& mesh, const int nCells, 
		     const VecD3Array& cellCentroid,
		     const string fileBase);

void markCell( Mesh& mesh, const int nCells, const int nSelfCells, 
	       const CRConnectivity& cellParticles, const CRConnectivity& cellCells );

const  shared_ptr<CRConnectivity> setParticleCells
                          (const StorageSite& rowSite,
			   const StorageSite& colSite, 
			   const Array<int> & connectivity);

const shared_ptr<CRConnectivity> setibFaceCells 
                          (const Mesh& mesh,
			   const Array<int>& ibFaceGroup,
			   const StorageSite& ibFaces, 
			   const StorageSite& cells,
			   const CRConnectivity& faceCells,
			   const CRConnectivity& cellFaces,
			   const VecD3Array& faceCentroid );

const shared_ptr<CRConnectivity> setibFaceParticles 
                          (const Mesh& mesh,
			   const StorageSite& ibFaces, 
			   const Array<int>& ibFaceGroup,
			   const StorageSite& particles,
			   const CRConnectivity& faceCells, 
			   const CRConnectivity& cellParticles,
			   const CRConnectivity& cellCells,
			   const Array<int>& particleTyp);

void markIBFaces(Mesh& mesh, const int nCells, 
		 const CRConnectivity& faceCells);

void checkIBFaces(const Array<int> & ibFaceList,
		  const VecD3Array& faceArea,
		  const CRConnectivity& faceCells,
		  const Mesh& mesh);

#endif
