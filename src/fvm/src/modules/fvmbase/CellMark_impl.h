// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CELLMARK_IMPL_H_
#define _CELLMARK_IMPL_H_

//this serves as a main function to mark the cells

#include "Field.h"
#include "Mesh.h"
#include "StorageSite.h"
#include "GeomFields.h"
#include "Vector.h"
#include "NumType.h"
#include "Octree.h"
#include "CellMark.h"
#include "Array.h"
#include "MPM_Particles.h"

 

void CellMark_Impl(Mesh& mesh, const GeomFields& geomFields, const string fileBase, 
		   Octree& O, MPM& solid, const int option);


 
#endif
