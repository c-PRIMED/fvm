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
class VacancyModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;
  typedef Array<Gradient<VectorT3> > VGradArray;
  
  Impl(const GeomFields& geomFields,
       VacancyFields& vacancyFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _vacancyFields(vacancyFields),
    _concentrationGradientModel(_meshes,_vacancyFields.concentration,
                              _vacancyFields.concentrationGradient,_geomFields),
    _fluxGradientModel(_meshes,
                       _vacancyFields.concentrationGradientVector,
                       _vacancyFields.plasticStrain,_geomFields),
    _initialNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        VacancyVC<T> *vc(new VacancyVC<T>());
        vc->vcType = "flow";
       _vcMap[mesh.getID()] = vc;
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            VacancyBC<T> *bc(new VacancyBC<T>());
            
            _bcMap[fg.id] = bc;

            if ((fg.groupType == "wall") ||
                (fg.groupType == "symmetry"))
            {
                bc->bcType = "SpecifiedVacaFlux";
            }
            else if ((fg.groupType == "velocity-inlet") ||
                     (fg.groupType == "pressure-outlet"))
            {
                bc->bcType = "SpecifiedConcentration";
            }
            else
              throw CException("VacancyModel: unknown face group type "
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
        const VacancyVC<T>& vc = *_vcMap[mesh.getID()];

	//concentration
        shared_ptr<TArray> tCell(new TArray(cells.getCountLevel1()));
        *tCell = _options["initialConcentration"];
        _vacancyFields.concentration.addArray(cells,tCell);
	
	if(_options.transient)
	  {
	    _vacancyFields.concentrationN1.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	    if (_options.timeDiscretizationOrder > 1)
	      _vacancyFields.concentrationN2.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	  }

	//diffusioncoefficient
        shared_ptr<TArray> condCell(new TArray(cells.getCountLevel1()));
        *condCell = vc["vacancyDiffusioncoefficient"];
        _vacancyFields.diffusioncoefficient.addArray(cells,condCell);
	
	//source 
	shared_ptr<TArray> sCell(new TArray(cells.getCountLevel1()));
	*sCell = T(0.);
	_vacancyFields.source.addArray(cells,sCell);

	//create a zero field
	shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
	*zeroCell = T(0.0);
	_vacancyFields.zero.addArray(cells,zeroCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCountLevel1()));
	*oneCell = T(1.0);
	_vacancyFields.one.addArray(cells,oneCell);

	//create specific vaca field   rho*Cp
	shared_ptr<TArray> cp(new TArray(cells.getCount()));
	*cp = vc["density"] * vc["specificVaca"];
	_vacancyFields.specificVaca.addArray(cells, cp);

        //initial temparature gradient array
        shared_ptr<TGradArray> gradT(new TGradArray(cells.getCountLevel1()));
        gradT->zero();
        _vacancyFields.concentrationGradient.addArray(cells,gradT);


        shared_ptr<VectorT3Array> vGrad(new VectorT3Array(cells.getCountLevel1()));
        vGrad->zero();
        _vacancyFields.concentrationGradientVector.addArray(cells,vGrad);



        shared_ptr<VGradArray> plasticStrainField(new VGradArray(cells.getCountLevel1()));
        plasticStrainField->zero();
        _vacancyFields.plasticStrain.addArray(cells,plasticStrainField);


	//inital convection flux at faces

	const StorageSite& allFaces = mesh.getFaces();
	shared_ptr<TArray> convFlux(new TArray(allFaces.getCount()));
	convFlux->zero();
	_vacancyFields.convectionFlux.addArray(allFaces,convFlux);

	//vaca flux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _vacancyFields.vacaFlux.addArray(faces,fluxFace);
          
        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _vacancyFields.vacaFlux.addArray(faces,fluxFace);
          
        }
	
	
    }
    _vacancyFields.diffusioncoefficient.syncLocal();
    _niters  =0;
    _initialNorm = MFRPtr();
  }
  
  VacancyBCMap& getBCMap() {return _bcMap;}
  VacancyVCMap& getVCMap() {return _vcMap;}

  VacancyBC<T>& getBC(const int id) {return *_bcMap[id];}

  VacancyModelOptions<T>& getOptions() {return _options;}

  void initLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_vacancyFields.concentration,&cells);

        ls.getX().addArray(tIndex,_vacancyFields.concentration.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_vacancyFields.vacaFlux,&faces);
            ls.getX().addArray(fIndex,_vacancyFields.vacaFlux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&_vacancyFields.vacaFlux,&faces);
            ls.getX().addArray(fIndex,_vacancyFields.vacaFlux.getArrayPtr(faces));

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
    _concentrationGradientModel.compute();
    
    
    DiscrList discretizations;
   
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _vacancyFields.concentration,
	  _vacancyFields.diffusioncoefficient,
	  _vacancyFields.concentrationGradient));
    discretizations.push_back(dd);
   
    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _vacancyFields.concentration,
	  _vacancyFields.convectionFlux,
	  _vacancyFields.zero,
	  _vacancyFields.concentrationGradient,
          _options.useCentralDifference));
    discretizations.push_back(cd);
    
    
    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  _vacancyFields.concentration,
	  _vacancyFields.source));
    discretizations.push_back(sd);
    
    if (_options.transient)
      {
	shared_ptr<Discretization>
	  td(new TimeDerivativeDiscretization<T, T, T>
	     (_meshes, _geomFields, 
	      _vacancyFields.concentration, 
	      _vacancyFields.concentrationN1,
	      _vacancyFields.concentrationN2,
	      _vacancyFields.specificVaca,
	      _options["timeStep"]));
	discretizations.push_back(td);
      }
    
    shared_ptr<Discretization>
      ibm(new GenericIBDiscretization<T,T,T>
	  (_meshes,_geomFields,_vacancyFields.concentration));
      
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

            const VacancyBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _vacancyFields.concentration,
                                  _vacancyFields.vacaFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedConcentration")
            {
                FloatValEvaluator<T>
                  bT(bc.getVal("specifiedConcentration"),faces);
                if (_vacancyFields.convectionFlux.hasArray(faces))
                {
                    const TArray& convectingFlux =
                      dynamic_cast<const TArray&>
                      (_vacancyFields.convectionFlux[faces]);
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
            else if (bc.bcType == "SpecifiedVacaFlux")
            {
                FloatValEvaluator<T>
                    bVacaFlux(bc.getVal("specifiedVacaFlux"),faces);
                    
                const int nFaces = faces.getCount();
                                
                for(int f=0; f<nFaces; f++)
                    {                        
                        gbc.applyNeumannBC(f, bVacaFlux[f]);
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
	        FloatValEvaluator<T> Xinf(bc.getVal("farFieldConcentration"), faces);
		const int nFaces = faces.getCount();
		for(int f=0; f<nFaces; f++)
		    gbc.applyConvectionBC(f, hCoeff[f], Xinf[f]);
	    }
            else
              throw CException(bc.bcType + " not implemented for VacancyModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _vacancyFields.concentration,
                                  _vacancyFields.vacaFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }
  
  T getVacaFluxIntegral(const Mesh& mesh, const int faceGroupId)
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
            const TArray& vacaFlux =
              dynamic_cast<const TArray&>(_vacancyFields.vacaFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += vacaFlux[f];
            found=true;
        }
    }
    if (!found)
      throw CException("getVacaFluxIntegral: invalid faceGroupID");
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
    foreach(typename VacancyBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename VacancyBC<T>::value_type& vp, *pos.second)
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
	
      TArray& concentration =
          dynamic_cast<TArray&>(_vacancyFields.concentration[cells]);
      TArray& concentrationN1 =
          dynamic_cast<TArray&>(_vacancyFields.concentrationN1[cells]);
     
      if (_options.timeDiscretizationOrder > 1)
        {
	  TArray& concentrationN2 =
	    dynamic_cast<TArray&>(_vacancyFields.concentrationN2[cells]);
	  concentrationN2 = concentrationN1;
        }
      concentrationN1 = concentration;
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
    
    MultiField::ArrayIndex tIndex(&_vacancyFields.concentration,&cells);

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
  
  void computeIBFaceConcentration(const StorageSite& particles)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;

    const TArray& pT =
      dynamic_cast<const TArray&>(_vacancyFields.concentration[particles]);

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
             dynamic_cast<const TArray&>(_vacancyFields.concentration[cells]);
    

           ibT->zero();

	   mIC.multiplyAndAdd(*ibT,cT);
	   mIP.multiplyAndAdd(*ibT,pT);
	  _vacancyFields.concentration.addArray(ibFaces,ibT);
        }  
    }

  }




  void computePlasticStrainRate()
  {
    const int numMeshes = _meshes.size();

    _concentrationGradientModel.compute();

    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      VectorT3Array& cGradVector = dynamic_cast<VectorT3Array&> (_vacancyFields.concentrationGradientVector[cells]);
      const TGradArray& cGrad = dynamic_cast<const TGradArray&> (_vacancyFields.concentrationGradient[cells]);
      for (int c=0; c<nCells; c++)
      {
        for (int k=0; k<3; k++)
          cGradVector[c][k] = cGrad[c][k];
      }
    }

    _fluxGradientModel.compute();


    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      VGradArray& plasticStrainCell = dynamic_cast<VGradArray&>(_vacancyFields.plasticStrain[cells]);
      const TArray& diffCoeff = dynamic_cast<const TArray&> (_vacancyFields.diffusioncoefficient[cells]);
      for(int c=0; c<nCells; c++)
      {
        Gradient<VectorT3>& ps = plasticStrainCell[c];
        
        const T df(diffCoeff[c]);
        Gradient<VectorT3> psTemp;
        for (int i=0;i<3;i++)
          for(int j=0;j<3;j++)
            psTemp[i][j] = 0.5 * df * (ps[i][j]+ps[j][i]);
        ps = psTemp;
      }
    }
  }





private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  VacancyFields& _vacancyFields;

  VacancyBCMap _bcMap;
  VacancyVCMap _vcMap;
  VacancyModelOptions<T> _options;
  GradientModel<T> _concentrationGradientModel;
  GradientModel<VectorT3> _fluxGradientModel;
  
  MFRPtr _initialNorm;
  int _niters;
};

template<class T>
VacancyModel<T>::VacancyModel(const GeomFields& geomFields,
                              VacancyFields& vacancyFields,
                              const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,vacancyFields,meshes))
{
  logCtor();
}


template<class T>
VacancyModel<T>::~VacancyModel()
{
  logDtor();
}

template<class T>
void
VacancyModel<T>::init()
{
  _impl->init();
}
 
template<class T>
void
VacancyModel<T>::computePlasticStrainRate()
{
  _impl->computePlasticStrainRate();
}

 
template<class T>
typename VacancyModel<T>::VacancyBCMap&
VacancyModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename VacancyModel<T>::VacancyVCMap&
VacancyModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
VacancyBC<T>&
VacancyModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
VacancyModelOptions<T>&
VacancyModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
VacancyModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
void
VacancyModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}

template<class T>
void
VacancyModel<T>::computeIBFaceConcentration(const StorageSite& particles)
{
  _impl->computeIBFaceConcentration(particles);
}

#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))

template<class T>
void
VacancyModel<T>::dumpMatrix(const string fileBase)
{
  _impl->dumpMatrix(fileBase);
}
#endif

template<class T>
T
VacancyModel<T>::getVacaFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
  return _impl->getVacaFluxIntegral(mesh, faceGroupId);
}

template<class T>
void
VacancyModel<T>::updateTime()
{
  _impl->updateTime();
}
