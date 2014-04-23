// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Mesh.h"
#include <sstream>

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
#include "ConvectionDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "SourceDiscretization.h"
#include "TimeDerivativeDiscretization.h"

template<class T>
class ThermalModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;
  
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
        ThermalVC<T> *vc(new ThermalVC<T>());
        vc->vcType = "flow";
       _vcMap[mesh.getID()] = vc;
        
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
        const ThermalVC<T>& vc = *_vcMap[mesh.getID()];

	//temperature
        shared_ptr<TArray> tCell(new TArray(cells.getCountLevel1()));
        *tCell = _options["initialTemperature"];
        _thermalFields.temperature.addArray(cells,tCell);
	
	if(_options.transient)
	  {
	    _thermalFields.temperatureN1.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	    if (_options.timeDiscretizationOrder > 1)
	      _thermalFields.temperatureN2.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	  }

	//conductivity
        shared_ptr<TArray> condCell(new TArray(cells.getCountLevel1()));
        *condCell = vc["thermalConductivity"];
        _thermalFields.conductivity.addArray(cells,condCell);
	
	//source 
	shared_ptr<TArray> sCell(new TArray(cells.getCountLevel1()));
	*sCell = T(0.);
	_thermalFields.source.addArray(cells,sCell);

	//create a zero field
	shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
	*zeroCell = T(0.0);
	_thermalFields.zero.addArray(cells,zeroCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCountLevel1()));
	*oneCell = T(1.0);
	_thermalFields.one.addArray(cells,oneCell);

	//create specific heat field   rho*Cp
	shared_ptr<TArray> cp(new TArray(cells.getCount()));
	*cp = vc["density"] * vc["specificHeat"];
	_thermalFields.specificHeat.addArray(cells, cp);

	//initial temparature gradient array
	shared_ptr<TGradArray> gradT(new TGradArray(cells.getCountLevel1()));
	gradT->zero();
	_thermalFields.temperatureGradient.addArray(cells,gradT);
        
	//inital convection flux at faces

	const StorageSite& allFaces = mesh.getFaces();
	shared_ptr<TArray> convFlux(new TArray(allFaces.getCount()));
	convFlux->zero();
	_thermalFields.convectionFlux.addArray(allFaces,convFlux);

	//heat flux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _thermalFields.heatFlux.addArray(faces,fluxFace);
          
        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _thermalFields.heatFlux.addArray(faces,fluxFace);
          
        }
	
	
    }
    _thermalFields.conductivity.syncLocal();
    _niters  =0;
    _initialNorm = MFRPtr();
  }
  
  ThermalBCMap& getBCMap() {return _bcMap;}
  ThermalVCMap& getVCMap() {return _vcMap;}

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

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
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
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _thermalFields.temperature,
	  _thermalFields.conductivity,
	  _thermalFields.temperatureGradient));
    discretizations.push_back(dd);
   
    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _thermalFields.temperature,
	  _thermalFields.convectionFlux,
	  _thermalFields.zero,
	  _thermalFields.temperatureGradient,
          _options.useCentralDifference));
    discretizations.push_back(cd);
    
    
    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  _thermalFields.temperature,
	  _thermalFields.source));
    discretizations.push_back(sd);
    
    if (_options.transient)
      {
	shared_ptr<Discretization>
	  td(new TimeDerivativeDiscretization<T, T, T>
	     (_meshes, _geomFields, 
	      _thermalFields.temperature, 
	      _thermalFields.temperatureN1,
	      _thermalFields.temperatureN2,
	      _thermalFields.specificHeat,
	      _options["timeStep"]));
	discretizations.push_back(td);
      }
    
    shared_ptr<Discretization>
      ibm(new GenericIBDiscretization<T,T,T>
	  (_meshes,_geomFields,_thermalFields.temperature));
      
    discretizations.push_back(ibm);
    

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
                FloatValEvaluator<T>
                  bT(bc.getVal("specifiedTemperature"),faces);
                if (_thermalFields.convectionFlux.hasArray(faces))
                {
                    const TArray& convectingFlux =
                      dynamic_cast<const TArray&>
                      (_thermalFields.convectionFlux[faces]);
                    const int nFaces = faces.getCount();
                                
                    for(int f=0; f<nFaces; f++)
                    {
                        if (convectingFlux[f] > 0.)
                        {
                            gbc.applyExtrapolationBC(f);
                        }
                        else
                        {
                            gbc.applyDirichletBC(f,bT[f]);
                        }
                    }
                }
                else
                    gbc.applyDirichletBC(bT);
            }
            else if (bc.bcType == "SpecifiedHeatFlux")
            {
                FloatValEvaluator<T>
                    bHeatFlux(bc.getVal("specifiedHeatFlux"),faces);
                    
                const int nFaces = faces.getCount();
                                
                for(int f=0; f<nFaces; f++)
                    {                        
                        gbc.applyNeumannBC(f, bHeatFlux[f]);
                    }                              
            }
            else if (bc.bcType == "Symmetry")
            {
                 T zeroFlux(NumTypeTraits<T>::getZero());
                 gbc.applyNeumannBC(zeroFlux);
            }
	    else if (bc.bcType == "Convective")
	    {
	        FloatValEvaluator<T> hCoeff(bc.getVal("convectiveCoefficient"), faces);
	        FloatValEvaluator<T> Xinf(bc.getVal("farFieldTemperature"), faces);
		const int nFaces = faces.getCount();
		for(int f=0; f<nFaces; f++)
		    gbc.applyConvectionBC(f, hCoeff[f], Xinf[f]);
	    }
	    else if (bc.bcType == "Radiative")
	    {
	        FloatValEvaluator<T> emissivity(bc.getVal("surfaceEmissivity"), faces);
	        FloatValEvaluator<T> Xinf(bc.getVal("farFieldTemperature"), faces);
		const int nFaces = faces.getCount();
		for(int f=0; f<nFaces; f++)
		    gbc.applyRadiationBC(f, emissivity[f], Xinf[f]);
	    }
	    else if (bc.bcType == "Mixed")
	    {
	        FloatValEvaluator<T> hCoeff(bc.getVal("convectiveCoefficient"), faces);
	        FloatValEvaluator<T> emissivity(bc.getVal("surfaceEmissivity"), faces);
	        FloatValEvaluator<T> Xinf(bc.getVal("farFieldTemperature"), faces);
	        const int nFaces = faces.getCount();
	        for(int f=0; f<nFaces; f++)
		    gbc.applyMixedBC(f, hCoeff[f], emissivity[f], Xinf[f]);
	    }
            else
              throw CException(bc.bcType + " not implemented for ThermalModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _thermalFields.temperature,
                                  _thermalFields.heatFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }
  
  T getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    T r(0.);
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const TArray& heatFlux =
              dynamic_cast<const TArray&>(_thermalFields.heatFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += heatFlux[f];
            found=true;
        }
    }
    if (!found)
      throw CException("getHeatFluxIntegral: invalid faceGroupID");
    return r;
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

        MFRPtr rNorm(_options.getLinearSolver().solve(ls));

        if (!_initialNorm) _initialNorm = rNorm;
        
        MFRPtr normRatio((*rNorm)/(*_initialNorm));

        cout << _niters << ": " << *rNorm << endl;

        
        _options.getLinearSolver().cleanup();

        ls.postSolve();
        ls.updateSolution();

        _niters++;
        if (*rNorm < _options.absoluteTolerance ||
            *normRatio < _options.relativeTolerance)
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
            cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
        }
    }
  }


  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    {
     
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
	
      TArray& temperature =
          dynamic_cast<TArray&>(_thermalFields.temperature[cells]);
      TArray& temperatureN1 =
          dynamic_cast<TArray&>(_thermalFields.temperatureN1[cells]);
     
      if (_options.timeDiscretizationOrder > 1)
        {
	  TArray& temperatureN2 =
	    dynamic_cast<TArray&>(_thermalFields.temperatureN2[cells]);
	  temperatureN2 = temperatureN1;
        }
      temperatureN1 = temperature;
    }
  }


#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
  
  void dumpMatrix(const string fileBase)
  {
    LinearSystem ls;
    initLinearization(ls);
    
    ls.initAssembly();
    
    linearize(ls);
    
    ls.initSolve();

    MultiFieldMatrix& matrix = ls.getMatrix();
    MultiField& b = ls.getB();
for ( unsigned int id = 0; id < _meshes.size(); id++ ){
    const Mesh& mesh = *_meshes[id];
    const StorageSite& cells = mesh.getCells();
    
    MultiField::ArrayIndex tIndex(&_thermalFields.temperature,&cells);

    T_Matrix& tMatrix =
      dynamic_cast<T_Matrix&>(matrix.getMatrix(tIndex,tIndex));

    TArray& tDiag = tMatrix.getDiag();
    TArray& tCoeff = tMatrix.getOffDiag();

    TArray& rCell = dynamic_cast<TArray&>(b[tIndex]);

    const CRConnectivity& cr = tMatrix.getConnectivity();

    const Array<int>& row = cr.getRow();
    const Array<int>& col = cr.getCol();
    
    const int nCells = cells.getSelfCount();
    int nCoeffs = nCells;

    for(int i=0; i<nCells; i++)
      for(int jp=row[i]; jp<row[i+1]; jp++)
      {
          const int j = col[jp];
          if (j<nCells) nCoeffs++;
      }
    stringstream ss;
    ss << id;
    string matFileName = fileBase + "_mesh" + ss.str() +  ".mat";


    FILE *matFile = fopen(matFileName.c_str(),"wb");
    
    fprintf(matFile,"%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matFile,"%d %d %d\n", nCells,nCells,nCoeffs);

    for(int i=0; i<nCells; i++)
    {
        fprintf(matFile,"%d %d %lf\n", i+1, i+1, tDiag[i]);
        for(int jp=row[i]; jp<row[i+1]; jp++)
        {
            const int j = col[jp];
            if (j<nCells)
              fprintf(matFile,"%d %d %lf\n", i+1, j+1, tCoeff[jp]);
        }
    }

    fclose(matFile);

    string rhsFileName = fileBase + ".rhs";
    FILE *rhsFile = fopen(rhsFileName.c_str(),"wb");
    
    for(int i=0; i<nCells; i++)
      fprintf(rhsFile,"%lf\n",-rCell[i]);

    fclose(rhsFile);

  }
}
#endif
  
  void computeIBFaceTemperature(const StorageSite& particles)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;

    const TArray& pT =
      dynamic_cast<const TArray&>(_thermalFields.temperature[particles]);

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
           const StorageSite& cells = mesh.getCells();
           const StorageSite& ibFaces = mesh.getIBFaces();
        
           GeomFields::SSPair key1(&ibFaces,&cells);
           const IMatrix& mIC =
           dynamic_cast<const IMatrix&>
           (*_geomFields._interpolationMatrices[key1]);

           GeomFields::SSPair key2(&ibFaces,&particles);
           const IMatrix& mIP =
           dynamic_cast<const IMatrix&>
           (*_geomFields._interpolationMatrices[key2]);

           shared_ptr<TArray> ibT(new TArray(ibFaces.getCount()));
        
           const TArray& cT =
             dynamic_cast<const TArray&>(_thermalFields.temperature[cells]);
    

           ibT->zero();

	   mIC.multiplyAndAdd(*ibT,cT);
	   mIP.multiplyAndAdd(*ibT,pT);
	  _thermalFields.temperature.addArray(ibFaces,ibT);
        }  
    }

  }

private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  ThermalFields& _thermalFields;

  ThermalBCMap _bcMap;
  ThermalVCMap _vcMap;
  ThermalModelOptions<T> _options;
  GradientModel<T> _temperatureGradientModel;
  
  MFRPtr _initialNorm;
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
typename ThermalModel<T>::ThermalVCMap&
ThermalModel<T>::getVCMap() {return _impl->getVCMap();}

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

template<class T>
void
ThermalModel<T>::computeIBFaceTemperature(const StorageSite& particles)
{
  _impl->computeIBFaceTemperature(particles);
}

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))

template<class T>
void
ThermalModel<T>::dumpMatrix(const string fileBase)
{
  _impl->dumpMatrix(fileBase);
}
#endif

template<class T>
T
ThermalModel<T>::getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
  return _impl->getHeatFluxIntegral(mesh, faceGroupId);
}

template<class T>
void
ThermalModel<T>::updateTime()
{
  _impl->updateTime();
}
