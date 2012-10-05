// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.


#include "Linearizer.h"
#include "Mesh.h"
#include "MultiFieldMatrix.h"
#include "MultiField.h"
#include "Discretization.h"
//#include <omp.h>

Linearizer::Linearizer()
{}

void
Linearizer::linearize(DiscrList& discretizations,
                      const MeshList& meshes, MultiFieldMatrix& matrix,
                      MultiField& x, MultiField& r)
{
  const int nMeshes = meshes.size();
  const int nDiscretizations = discretizations.size();
  
  //#pragma omp parallel for
  for(int n=0; n<nMeshes; n++)
  {
      const Mesh& mesh = *meshes[n];
      if (mesh.isShell() == false){
	for(int nd=0; nd<nDiscretizations; nd++)
	  discretizations[nd]->discretize(mesh,matrix,x,r);
      }
  }
}

