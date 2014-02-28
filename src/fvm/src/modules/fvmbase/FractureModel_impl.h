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
//#include "ConvectionDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
//#include "GenericIBDiscretization.h"
#include "SourceDiscretizationforFracture.h"
#include "TimeDerivativeDiscretization.h"
#include "SquareTensor.h"

template<class T>
class FractureModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;
  typedef SquareTensor<T,3>  DiagTensorT3;
  
  Impl(const GeomFields& geomFields,
       FractureFields& fractureFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _fractureFields(fractureFields),
    _phasefieldGradientModel(_meshes,_fractureFields.phasefieldvalue,
                              _fractureFields.phasefieldGradient,_geomFields),
    _initialNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        FractureVC<T> *vc(new FractureVC<T>());
        vc->vcType = "flow";
       _vcMap[mesh.getID()] = vc;
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            FractureBC<T> *bc(new FractureBC<T>());
            
            _bcMap[fg.id] = bc;

            if ((fg.groupType == "wall") ||
                (fg.groupType == "symmetry"))
            {
                bc->bcType = "SpecifiedPhaseFieldFlux";
            }
            else if ((fg.groupType == "velocity-inlet") ||
                     (fg.groupType == "pressure-outlet"))
            {
                bc->bcType = "SpecifiedPhaseFieldValue";
            }
            else
              throw CException("FractureModel: unknown face group type "
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
        const FractureVC<T>& vc = *_vcMap[mesh.getID()];

	//phasefieldvalue
        shared_ptr<TArray> tCell(new TArray(cells.getCountLevel1()));
        *tCell = _options["initialPhaseFieldValue"];
        _fractureFields.phasefieldvalue.addArray(cells,tCell);
	
	if(_options.transient)
	  {
	    _fractureFields.phasefieldvalueN1.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	    if (_options.timeDiscretizationOrder > 1)
	      _fractureFields.phasefieldvalueN2.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	  }

	//conductivity
        shared_ptr<TArray> condCell(new TArray(cells.getCountLevel1()));
        *condCell = vc["fractureConductivity"];
        _fractureFields.conductivity.addArray(cells,condCell);
	
	//source 
	shared_ptr<TArray> sCell(new TArray(cells.getCountLevel1()));
	*sCell = vc["fractureSource"];
	//*sCell =T(0.0);
	_fractureFields.source.addArray(cells,sCell);

	//source coef
	shared_ptr<TArray> scoefCell(new TArray(cells.getCountLevel1()));
	*scoefCell = vc["fractureSourceCoef"];
	//*scoefCell=T(0.0);
	_fractureFields.sourcecoef.addArray(cells,scoefCell);

	//create a zero field
	shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
	*zeroCell = T(0.0);
	_fractureFields.zero.addArray(cells,zeroCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCountLevel1()));
	*oneCell = T(1.0);
	_fractureFields.one.addArray(cells,oneCell);

	//create specific heat field   rho*Cp
	//shared_ptr<TArray> cp(new TArray(cells.getCount()));
	//*cp = vc["density"] * vc["specificHeat"];
	//_thermalFields.specificHeat.addArray(cells, cp);

	//initial phasefieldvalue gradient array
	shared_ptr<TGradArray> gradT(new TGradArray(cells.getCountLevel1()));
	gradT->zero();
	_fractureFields.phasefieldGradient.addArray(cells,gradT);
        
	//inital convection flux at faces
	
	//const StorageSite& allFaces = mesh.getFaces();
	//shared_ptr<TArray> convFlux(new TArray(allFaces.getCount()));
	//convFlux->zero();
	//_thermalFields.convectionFlux.addArray(allFaces,convFlux);

	//phasefield flux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _fractureFields.phasefieldFlux.addArray(faces,fluxFace);
          
        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _fractureFields.phasefieldFlux.addArray(faces,fluxFace);
          
        }
	
	
    }
    _fractureFields.conductivity.syncLocal();
    _niters  =0;
    _initialNorm = MFRPtr();
  }
  
  FractureBCMap& getBCMap() {return _bcMap;}
  FractureVCMap& getVCMap() {return _vcMap;}

  FractureBC<T>& getBC(const int id) {return *_bcMap[id];}

  FractureModelOptions<T>& getOptions() {return _options;}

  void initLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_fractureFields.phasefieldvalue,&cells);

        ls.getX().addArray(tIndex,_fractureFields.phasefieldvalue.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_fractureFields.phasefieldFlux,&faces);
            ls.getX().addArray(fIndex,_fractureFields.phasefieldFlux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&_fractureFields.phasefieldFlux,&faces);
            ls.getX().addArray(fIndex,_fractureFields.phasefieldFlux.getArrayPtr(faces));

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
    _phasefieldGradientModel.compute();
    
    DiscrList discretizations;

    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _fractureFields.phasefieldvalue,
	  _fractureFields.conductivity,
	  _fractureFields.phasefieldGradient));
    discretizations.push_back(dd);
 
    //shared_ptr<Discretization>
    //  cd(new ConvectionDiscretization<T,T,T>
	//(_meshes,_geomFields,
	//  _fractureFields.phasefieldvalue,
	//  _fractureFields.phasefieldFlux,
	//  _fractureFields.zero,
	//  _fractureFields.phasefieldGradient,
    //      _options.useCentralDifference));
    //discretizations.push_back(cd);
    
    shared_ptr<Discretization>
      scfd(new SourceDiscretizationforFracture<T, T, T>
	 (_meshes, 
	  _geomFields, 
	  _fractureFields.phasefieldvalue,
	  _fractureFields.source,
	  _fractureFields.sourcecoef));
    discretizations.push_back(scfd);
    
    if (_options.transient)
      {
	shared_ptr<Discretization>
	  td(new TimeDerivativeDiscretization<T, T, T>
	     (_meshes, _geomFields, 
	      _fractureFields.phasefieldvalue, 
	      _fractureFields.phasefieldvalueN1,
	      _fractureFields.phasefieldvalueN2,
	      _fractureFields.one,
	      _options["timeStep"]));
	discretizations.push_back(td);
      }
    
    //shared_ptr<Discretization>
    //  ibm(new GenericIBDiscretization<T,T,T>
	//  (_meshes,_geomFields,_thermalFields.temperature));
    //discretizations.push_back(ibm);
    

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

            const FractureBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _fractureFields.phasefieldvalue,
                                  _fractureFields.phasefieldFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedPhaseFieldValue")
            {
	        FloatValEvaluator<T>
                  bT(bc.getVal("specifiedPhaseFieldValue"),faces);
		gbc.applyDirichletBC(bT);
            }
            else if (bc.bcType == "SpecifiedPhaseFieldFlux")
            {
                FloatValEvaluator<T>
                    bPhaseFieldFlux(bc.getVal("specifiedPhaseFieldFlux"),faces);
                    
                const int nFaces = faces.getCount();
                                
                for(int f=0; f<nFaces; f++)
                    {                        
                        gbc.applyNeumannBC(f, bPhaseFieldFlux[f]);
                    }                              
            }
            else if (bc.bcType == "Symmetry")
            {
                 T zeroFlux(NumTypeTraits<T>::getZero());
                 gbc.applyNeumannBC(zeroFlux);
            }
	    else
              throw CException(bc.bcType + " not implemented for FractureModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _fractureFields.phasefieldvalue,
                                  _fractureFields.phasefieldFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }
  
 /* T getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
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
  }	*/
	

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
    foreach(typename FractureBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename FractureBC<T>::value_type& vp, *pos.second)
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
	
      TArray& phasefieldvalue =
          dynamic_cast<TArray&>(_fractureFields.phasefieldvalue[cells]);
      TArray& phasefieldvalueN1 =
          dynamic_cast<TArray&>(_fractureFields.phasefieldvalueN1[cells]);
     
      if (_options.timeDiscretizationOrder > 1)
        {
	  TArray& phasefieldvalueN2 =
	    dynamic_cast<TArray&>(_fractureFields.phasefieldvalueN2[cells]);
	  phasefieldvalueN2 = phasefieldvalueN1;
        }
      phasefieldvalueN1 = phasefieldvalue;
    }
  }

/*
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

  }	*/

private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  FractureFields& _fractureFields;

  FractureBCMap _bcMap;
  FractureVCMap _vcMap;
  FractureModelOptions<T> _options;
  GradientModel<T> _phasefieldGradientModel;
  
  MFRPtr _initialNorm;
  int _niters;
};

template<class T>
FractureModel<T>::FractureModel(const GeomFields& geomFields,
                              FractureFields& fractureFields,
                              const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,fractureFields,meshes))
{
  logCtor();
}


template<class T>
FractureModel<T>::~FractureModel()
{
  logDtor();
}

template<class T>
void
FractureModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename FractureModel<T>::FractureBCMap&
FractureModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename FractureModel<T>::FractureVCMap&
FractureModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
FractureBC<T>&
FractureModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
FractureModelOptions<T>&
FractureModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
FractureModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
void
FractureModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}

/*
template<class T>
void
FractureModel<T>::computeIBFaceTemperature(const StorageSite& particles)
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
}	*/

template<class T>
void
FractureModel<T>::updateTime()
{
  _impl->updateTime();
}
