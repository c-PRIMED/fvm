// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Mesh.h"
#include "ThermalFields.h"
#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "StorageSite.h"

template<class T>
class RosselandModel<T>::Impl
{
public:
  typedef Array<T> TArray;

  Impl(const GeomFields& geomFields,
       ThermalFields& thermalFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _thermalFields(thermalFields)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        RosselandVC<T> *vc(new RosselandVC<T>());
        _vcMap[mesh.getID()] = vc;
    }
  }

  VCMap& getVCMap() {return _vcMap;}

  void init() { advance(1,true); }

  bool advance(const int niter, bool init = false)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const RosselandVC<T>& vc = *_vcMap[mesh.getID()];

        const StorageSite& cells = mesh.getCells();

        TArray& conductivity= dynamic_cast<TArray&>(_thermalFields.conductivity[cells]);

        FloatValEvaluator<T>  cellT(vc.getVal("temperature"),cells);

        const int nCells = cells.getCount();
        const float stefan = 0.0000000567;
        const T index(1.); 
        const T a(1.);
        const T sigmas(0.);
        const T C(1.0);
        const T gamma = 1/(3*(a+sigmas) - C*sigmas);
        for(int c=0; c<nCells; c++)
        {
            
            T temp = cellT[c];
                        
	    conductivity[c] = 16*stefan*gamma*index*index*temp*temp*temp;
        }
    }
    return true;
  }


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  ThermalFields& _thermalFields;

  VCMap _vcMap;
};

template<class T>
RosselandModel<T>::RosselandModel(const GeomFields& geomFields,
                        ThermalFields& thermalFields,
                        const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,thermalFields,meshes))
{
  logCtor();
}


template<class T>
RosselandModel<T>::~RosselandModel()
{
  logDtor();
}

template<class T>
void
RosselandModel<T>::init()
{
  _impl->init();
}

template<class T>
typename RosselandModel<T>::VCMap&
RosselandModel<T>::getVCMap() {return _impl->getVCMap();}


template<class T>
bool
RosselandModel<T>::advance(const int niter)
{
  return _impl->advance(niter);
}

   

