// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include "LinearSystem.h"
#include "Array.h"
#include "Field.h"

LinearSystem::LinearSystem() :
  isSymmetric(false),
  _x(new MultiField()),
  _b(),
  _delta(),
  _residual(),
  _coarseIndex(),
  _coarseningField(0)
{}

LinearSystem::~LinearSystem()
{
  

}



void
LinearSystem::initAssembly()
{
  _matrix.initAssembly();
  _b = dynamic_pointer_cast<MultiField>(_x->newClone());
  _b->zero();
}

void
LinearSystem::initSolve()
{
  _delta = dynamic_pointer_cast<MultiField>(_x->newClone());
  _residual = dynamic_pointer_cast<MultiField>(_x->newClone());
  _delta->zero();
  _residual->zero();


  MultiField::ArrayIndexList indicesToRemove;

  const MultiField::ArrayIndexList& arrayIndices = _b->getArrayIndices();
  foreach(MultiField::ArrayIndex i, arrayIndices)
  {
      if (!_matrix.hasMatrix(i,i) || _matrix.getMatrix(i,i).isInvertible())
        indicesToRemove.push_back(i);
      else
        _matrix.getMatrix(i,i).eliminateBoundaryEquations((*_b)[i]);
  }

  _xAux = _x->extract(indicesToRemove);
  _bAux = _b->extract(indicesToRemove);
  _deltaAux = _delta->extract(indicesToRemove);
  _residualAux = _residual->extract(indicesToRemove);
  
}

shared_ptr<LinearSystem>
LinearSystem::createCoarse(const int groupSize, const double weightRatioThreshold)
{

  shared_ptr<LinearSystem> coarseLS(new LinearSystem());
  /**
   * we create only one entry in coarseIndex for each
   * StorageSite even if the site is present in multiple fine
   * ArrayIndices of x. If a coarseningField is specified we will
   * use that otherwise the first one gets used
   * 
   */
 
  // keep track of which fieldIndex was used for coarsening each site
  map<const StorageSite*, const Field *> sitesCoarsenedWithField; 
  

  const MultiField::ArrayIndexList& arrayIndices = _b->getArrayIndices();
  foreach(MultiField::ArrayIndex ai,arrayIndices)
  {
      const Field *fieldIndex = ai.first;
      const StorageSite *site = ai.second;

      if (!_coarseningField || fieldIndex == _coarseningField)
      {
          if (sitesCoarsenedWithField.find(site) == sitesCoarsenedWithField.end())
          {
              sitesCoarsenedWithField[site] = fieldIndex;
              const ArrayBase& bi = (*_b)[ai];
              shared_ptr<Array<int> > cIndex(new Array<int>(bi.getLength()));
              _coarseIndex.addArray(ai,cIndex);

          }
      }

  }
  
  /**
   *  coarsen the matrices corresponding to the fieldIndices we
   * have selected and then sync ghost cell values. This
   * process will fill up the coarseSizes and coarseGhostSizes
   * variables in the matrix as well create the mappers
   * 
   */

  _matrix.createCoarsening(_coarseIndex,groupSize,weightRatioThreshold);

  _coarseIndex.sync();


_matrix.syncGhostCoarsening(_coarseIndex);


  // we can now create the coarse sites for each fine site
  foreach(MultiField::ArrayIndex k,arrayIndices)
  {
      // we will only have entries for indices which were selected to be coarsened
      if (_matrix._coarseSizes.find(k) != _matrix._coarseSizes.end())
      {
          int selfSize = _matrix._coarseSizes[k];
          int ghostSize = _matrix._coarseGhostSizes[k];
          _matrix._coarseSites[k] =
            shared_ptr<StorageSite>(new StorageSite(selfSize,ghostSize));
      }
  }

  // set the mappers in the newly created coarse sites
  foreach(MultiField::ArrayIndex k,arrayIndices)
  {
      if (_matrix._coarseSites.find(k) != _matrix._coarseSites.end())
      {
          const StorageSite& fineSite = *k.second;

          StorageSite& coarseSite = *_matrix._coarseSites[k];

          const StorageSite::ScatterMap& fineScatterMap = fineSite.getScatterMap();
          StorageSite::ScatterMap& coarseScatterMap = coarseSite.getScatterMap();

          foreach(const StorageSite::ScatterMap::value_type& pos, fineScatterMap)
          {
              const StorageSite& fineOSite = *pos.first;
              MultiField::ArrayIndex kO(k.first,&fineOSite);

  
#ifdef FVM_PARALLEL
              // the ghost site will not have its corresponding coarse
              // site created yet so we create it here
              if (_matrix._coarseSites.find(kO) == _matrix._coarseSites.end())
              {
                  shared_ptr<StorageSite> ghostSite
                    (new StorageSite(-1));
                  ghostSite->setGatherProcID ( fineOSite.getGatherProcID() );
                  ghostSite->setScatterProcID( fineOSite.getScatterProcID() );
                  ghostSite->setTag( fineOSite.getTag() );
                  _matrix._coarseSites[kO] = ghostSite;
              }
#endif
              const StorageSite& coarseOSite = *_matrix._coarseSites[kO];
              MultiFieldMatrix::SSPair sskey(&fineSite,&fineOSite);
              coarseScatterMap[&coarseOSite] = _matrix._coarseScatterMaps[sskey];
          }

          const StorageSite::GatherMap& fineGatherMap = fineSite.getGatherMap();
          StorageSite::GatherMap& coarseGatherMap = coarseSite.getGatherMap();
          foreach(const StorageSite::GatherMap::value_type& pos, fineGatherMap)
          {
              const StorageSite& fineOSite = *pos.first;
              MultiField::ArrayIndex kO(k.first,&fineOSite);
              const StorageSite& coarseOSite = *_matrix._coarseSites[kO];
              MultiFieldMatrix::SSPair sskey(&fineSite,&fineOSite);


              coarseGatherMap[&coarseOSite] = _matrix._coarseGatherMaps[sskey];
          }

 	
      }
  }

  // create the connectivities for the coarse matrices
  _matrix.createCoarseToFineMapping(_coarseIndex);

  /**
   *  now we have the coarse indices, sizes, sites, mappers and *
   *  connectivities computed for the fieldIndices we *
   *  selected. Before we can create coarse matrices we need to * set
   *  these for the other fieldIndices since the matrices of * course
   *  need to be created for all the fine matrices * present
   * 
   */

  foreach(MultiField::ArrayIndex k,arrayIndices)
  {
      const StorageSite *site = k.second;
      if (_matrix._coarseSites.find(k) == _matrix._coarseSites.end())
      {
          const Field* fieldIndexUsedForCoarsening = sitesCoarsenedWithField[site];
          MultiField::ArrayIndex indexCoarsened(fieldIndexUsedForCoarsening,site);
          _coarseIndex.addArray(k, _coarseIndex.getArrayPtr(indexCoarsened));
          _matrix._coarseSizes[k] = _matrix._coarseSizes[indexCoarsened];
          _matrix._coarseGhostSizes[k] = _matrix._coarseGhostSizes[indexCoarsened];
          _matrix._coarseSites[k] = _matrix._coarseSites[indexCoarsened];
          _matrix._coarseToFineMappings[k] = _matrix._coarseToFineMappings[indexCoarsened];
      }
  }

  _matrix.createCoarseConnectivity(_coarseIndex);
  _matrix.createCoarseMatrices(_coarseIndex);
  foreach(MultiField::ArrayIndex fineRowIndex,arrayIndices)
  {
      MultiField::ArrayIndex coarseRowIndex (fineRowIndex.first,
                                             _matrix._coarseSites[fineRowIndex].get());
      foreach(MultiField::ArrayIndex fineColIndex,arrayIndices)
      {
          MultiFieldMatrix::EntryIndex fineEntryIndex(fineRowIndex,fineColIndex);
          if (_matrix.hasMatrix(fineRowIndex,fineColIndex))
          {
              MultiField::ArrayIndex coarseColIndex(fineColIndex.first,
                                                    _matrix._coarseSites[fineColIndex].get());
              MultiFieldMatrix::EntryIndex coarseEntryIndex(coarseRowIndex,coarseColIndex);
              coarseLS->_matrix._matrices[coarseEntryIndex] = 
                _matrix._coarseMatrices[fineEntryIndex];
          }
      }
  }

  coarseLS->_b = shared_ptr<MultiField>(new MultiField());

  foreach(MultiField::ArrayIndex k,arrayIndices)
  {
      const StorageSite& coarseSite = *_matrix._coarseSites[k];
      MultiField::ArrayIndex coarseIndex(k.first,&coarseSite);
      coarseLS->_b->addArray(coarseIndex, (*_b)[k].newSizedClone(coarseSite.getCount()));
  }

  coarseLS->_b->zero();
  coarseLS->_delta = dynamic_pointer_cast<MultiField>(coarseLS->_b->newClone());
  coarseLS->_delta->zero();
  coarseLS->_residual = dynamic_pointer_cast<MultiField>(coarseLS->_b->newClone());
  coarseLS->_residual->zero();

  return coarseLS;
}

void LinearSystem::postSolve()
{
  _matrix.solveBoundary(*_delta,*_b,*_residual);

  if (_xAux)
  {
      _x->merge(*_xAux);
      _b->merge(*_bAux);
      _delta->merge(*_deltaAux);
      _residual->merge(*_residualAux);

      _matrix.solveBoundary(*_delta,*_b,*_residual);
  }
}

void LinearSystem::updateSolution()
{
  *_x += *_delta;
   _x->sync();
}

