#include "Mesh.h"

#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
//#include "FieldSet.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"
#include "GenericBCS.h"
#include "Vector.h"
#include "DiffusionDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"

template<class T>
class ThermalModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  Impl(const GeomFields& geomFields,
       ThermalFields& thermalFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _thermalFields(thermalFields),
    _temperatureGradientModel(_meshes,_thermalFields.temperature,
                              _thermalFields.temperatureGradient,_geomFields),
    _initialNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            ThermalBC<T> *bc(new ThermalBC<T>());
            
            _bcMap[fg.id] = bc;
            if ((fg.groupType == "wall") ||
                (fg.groupType == "symmetry"))
            {
                bc->bcType = "SpecifiedHeatFlux";
            }
            else if ((fg.groupType == "velocity-inlet") ||
                     (fg.groupType == "pressure-outlet"))
            {
                bc->bcType = "SpecifiedTemperature";
            }
            else
              throw CException("ThermalModel: unknown face group type "
                               + fg.groupType);
        }
    }
  }

  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();

        shared_ptr<TArray> tCell(new TArray(cells.getCount()));

        *tCell = _options["initialTemperature"];
      
        _thermalFields.temperature.addArray(cells,tCell);
      
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _thermalFields.heatFlux.addArray(faces,fluxFace);
          
        }
    }

    _niters  =0;
    _initialNorm = MFPtr();
  }
  
  ThermalBCMap& getBCMap() {return _bcMap;}

  ThermalBC<T>& getBC(const int id) {return *_bcMap[id];}

  ThermalModelOptions<T>& getOptions() {return _options;}

  void initLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

        ls.getX().addArray(tIndex,_thermalFields.temperature.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_thermalFields.heatFlux,&faces);
            ls.getX().addArray(fIndex,_thermalFields.heatFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
    }
  }

  void linearize(LinearSystem& ls)
  {
    _temperatureGradientModel.compute();
    
    DiscrList discretizations;
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>(_meshes,_geomFields,
                                            _thermalFields.temperature,
                                            _thermalFields.conductivity,
                                            _thermalFields.temperatureGradient));
    discretizations.push_back(dd);
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            const ThermalBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _thermalFields.temperature,
                                  _thermalFields.heatFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedTemperature")
            {
                const T bT(bc["specifiedTemperature"]);
                gbc.applyDirichletBC(bT);
            }
            else if (bc.bcType == "SpecifiedHeatFlux")
            {
                const T specifiedFlux(bc["specifiedHeatFlux"]);
                gbc.applyNeumannBC(specifiedFlux);
            }
            else
              throw CException(bc.bcType + " not implemented for ThermalModel");
        }
    }
  }
  

  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)
    { 
        LinearSystem ls;
        initLinearization(ls);
        
        ls.initAssembly();

        linearize(ls);

        ls.initSolve();

        MFPtr rNorm(ls.getResidual().getOneNorm());

        if (!_initialNorm) _initialNorm = rNorm;
        
        MFPtr normRatio((*rNorm)/(*_initialNorm));

        cout << _niters << ": " << *rNorm << endl;


        AMG solver(ls);
        solver.solve();
        //        solver.cleanUp();

        ls.postSolve();
        ls.updateSolution();

        _niters++;
        if (*rNorm < _options["absoluteTolerance"] ||
            *normRatio < _options["relativeTolerance"])
          break;
    }
  }
    
  void printBCs()
  {
    foreach(typename ThermalBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename ThermalBC<T>::value_type& vp, *pos.second)
        {
            cout << "   " << vp.first << " "  << vp.second <<  endl;
        }
    }
  }
private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  ThermalFields& _thermalFields;

  ThermalBCMap _bcMap;
  ThermalModelOptions<T> _options;
  GradientModel<T> _temperatureGradientModel;
  
  MFPtr _initialNorm;
  int _niters;
};

template<class T>
ThermalModel<T>::ThermalModel(const GeomFields& geomFields,
                              ThermalFields& thermalFields,
                              const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,thermalFields,meshes))
{
  logCtor();
}


template<class T>
ThermalModel<T>::~ThermalModel()
{
  logDtor();
}

template<class T>
void
ThermalModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename ThermalModel<T>::ThermalBCMap&
ThermalModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
ThermalBC<T>&
ThermalModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
ThermalModelOptions<T>&
ThermalModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
ThermalModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
void
ThermalModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}
