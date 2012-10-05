// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <iostream>
#include "FVMParticles.h"
#include "CRConnectivity.h"

using namespace std;

FVMParticles::FVMParticles( const MeshList& meshList )
:_meshList( meshList ), _sweepIter(5)
{
   _nmesh     = _meshList.size();
   _cellIDSet.resize( _nmesh );

}


FVMParticles::~FVMParticles()
{
   
}

void   
FVMParticles::setParticles(int nsweep)
{
     _sweepIter = nsweep;
      assert( _sweepIter > 0 );
     //loop over meshes
     for ( int id = 0; id < _nmesh; id++){
        const Array<int>&  cellType   = _meshList.at(id)->getIBType();
        const Array<int>&  faceIBList = _meshList.at(id)->getIBFaceList();
        const CRConnectivity& cellCells = _meshList.at(id)->getCellCells();
        const CRConnectivity& faceCells  = _meshList.at(id)->getAllFaceCells();
        int num_inner_cells = _meshList.at(id)->getCells().getSelfCount();
        vector<int> sweep_particles_old;
        vector<int> sweep_particles_new;
        for ( int sweep = 0; sweep < _sweepIter; sweep++ ){
           //loop over cells on a mesh
           if ( sweep == 0 ) 
               for ( int n = 0; n < faceIBList.getLength(); n++){
                   //loop over surrounding cells of a cell only if it is immersed boundary or if it is a fvm particle
                   int  face_id = faceIBList[n];
                   int  cell0 = faceCells( face_id, 0);
                   int  cell1 = faceCells( face_id, 1);
                   int cell_id = -1;
                   if  ( cellType[cell0] == Mesh::IBTYPE_BOUNDARY  )  
                        cell_id = cell0;
                   else 
                        cell_id = cell1;

                   //cout << " immersed_cells id =  " << cell_id << endl;
                   int count_neigh_cells = cellCells.getCount(cell_id);
                   for ( int j = 0; j < count_neigh_cells; j++){
                       int neigh_cell_id = cellCells(cell_id, j);
                       bool is_fvm_particle = (_cellIDSet.at(id).count(neigh_cell_id) == 0) &&  (cellType[neigh_cell_id] == Mesh::IBTYPE_FLUID) && 
                                               (neigh_cell_id < num_inner_cells);
                       //accept only if it is a fluid and not in _cellIDSet
                       if ( is_fvm_particle ){
                          _cellIDSet.at(id).insert( neigh_cell_id );
                          sweep_particles_old.push_back( neigh_cell_id );
                       }
                    }  
              }

          if ( sweep > 0 )
               for ( int n = 0; n < int( sweep_particles_old.size() ); n++){
                   int cell_id = sweep_particles_old.at(n);
                   int count_neigh_cells = cellCells.getCount(n);
                   for ( int j = 0; j < count_neigh_cells; j++){
                       int neigh_cell_id = cellCells(cell_id, j);
                       bool is_fvm_particle =  (_cellIDSet.at(id).count(neigh_cell_id) == 0) &&  (cellType[neigh_cell_id] == Mesh::IBTYPE_FLUID) &&
                                              (neigh_cell_id < num_inner_cells);
                       //accept only if it is a fluid and not in _cellIDSet
                       if ( is_fvm_particle ){
                          _cellIDSet.at(id).insert( neigh_cell_id );
                          sweep_particles_new.push_back( neigh_cell_id );
                       }
                    }
               }

          //replace old with new one
          if ( sweep > 0 ){
             sweep_particles_old.resize( sweep_particles_new.size() );
             sweep_particles_old = sweep_particles_new;
          }
         cout << "sweep_end = " << sweep << endl;
        }
    }

    //setting array values from data sets
    for ( int id = 0; id < _nmesh; id++ ){
       int array_size = _cellIDSet.at(id).size();
      _cellID.push_back( ArrayIntPtr( new Array<int>(array_size) ) );
      set<int>::const_iterator it;
      int indx = 0;
      for ( it = _cellIDSet.at(id).begin(); it != _cellIDSet.at(id).end(); it++)
          (*_cellID.at(id))[indx++] = *it;
   } 


}

