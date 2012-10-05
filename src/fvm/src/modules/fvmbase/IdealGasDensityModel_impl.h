// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "StorageSite.h"

template<class T>
class IdealGasDensityModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  
  Impl(const GeomFields& geomFields,
       FlowFields& flowFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _flowFields(flowFields)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        IdealGasVC<T> *vc(new IdealGasVC<T>());
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

        const IdealGasVC<T>& vc = *_vcMap[mesh.getID()];
            
        const StorageSite& cells = mesh.getCells();

        TArray& density = dynamic_cast<TArray&>(_flowFields.density[cells]);
    
        FloatValEvaluator<T>  cellP(vc.getVal("pressure"),cells);
        FloatValEvaluator<T>  cellT(vc.getVal("temperature"),cells);

        const T operatingPressure = vc["operatingPressure"];
        const T molWt = vc["molecularWeight"];

        const T Rgas = 8314.472/molWt;
        const int nCells = cells.getCount();

        const T pMin(1000);
        const T TMin(1);
        const T urf = init ? T(1) : T(vc["urf"]);
        
        for(int c=0; c<nCells; c++)
        {
            T absP = (cellP[c] + operatingPressure);
            T temp = cellT[c];
            if (absP < pMin) absP  = pMin;
            if (temp < TMin) temp  = TMin;

            if (init)
              density[c] = absP/(Rgas*temp);
            else
              density[c] = urf*absP/(Rgas*temp) + (1.0-urf)*density[c];
        }
    }
    return true;
  }

  
private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  FlowFields& _flowFields;

  VCMap _vcMap;
};

template<class T>
IdealGasDensityModel<T>::IdealGasDensityModel(const GeomFields& geomFields,
                        FlowFields& flowFields,
                        const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,flowFields,meshes))
{
  logCtor();
}


template<class T>
IdealGasDensityModel<T>::~IdealGasDensityModel()
{
  logDtor();
}

template<class T>
void
IdealGasDensityModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename IdealGasDensityModel<T>::VCMap&
IdealGasDensityModel<T>::getVCMap() {return _impl->getVCMap();}


template<class T>
bool
IdealGasDensityModel<T>::advance(const int niter)
{
  return _impl->advance(niter);
}

