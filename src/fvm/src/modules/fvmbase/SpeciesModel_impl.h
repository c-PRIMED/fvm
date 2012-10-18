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
#include "TimeDerivativeDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "SourceDiscretization.h"
#include "LinearizeInterfaceJump.h"
#include "LinearizeSpeciesInterface.h"

template<class T>
class SpeciesModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;
  
  Impl(const GeomFields& geomFields,
       const MeshList& meshes,
       const int nSpecies) :
    _meshes(meshes),
    _geomFields(geomFields),
    _niters(0),
    _nSpecies(nSpecies),
    _speciesModelFields("speciesModel")
  {
    
    const int numMeshes = _meshes.size();

    for (int m=0; m<_nSpecies; m++)
    {
      SpeciesVCMap *svcmap = new SpeciesVCMap();
      _vcMapVector.push_back(svcmap);
      SpeciesBCMap *sbcmap = new SpeciesBCMap();
      _bcMapVector.push_back(sbcmap);
      SpeciesFields *sFields = new SpeciesFields("species");
      _speciesFieldsVector.push_back(sFields);
      MFRPtr *iNorm = new MFRPtr();
      _initialNormVector.push_back(iNorm);
      MFRPtr *rCurrent = new MFRPtr();
      _currentResidual.push_back(rCurrent); 
        
      for (int n=0; n<numMeshes; n++)
      {

        const Mesh& mesh = *_meshes[n];
        SpeciesVC<T> *vc(new SpeciesVC<T>());
        vc->vcType = "flow";
        (*svcmap)[mesh.getID()] = vc;
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            SpeciesBC<T> *bc(new SpeciesBC<T>());
            
            (*sbcmap)[fg.id] = bc;

            if ((fg.groupType == "wall") ||
                (fg.groupType == "symmetry"))
            {
                bc->bcType = "SpecifiedMassFlux";
            }
            else if ((fg.groupType == "velocity-inlet") ||
                     (fg.groupType == "pressure-outlet"))
            {
                bc->bcType = "SpecifiedMassFraction";
            }
            else
              throw CException("SpeciesModel: unknown face group type "
                               + fg.groupType);
        }
      }
    }
  }

  void init()
  {
    const int numMeshes = _meshes.size();

    for (int m=0; m<_nSpecies; m++)
    {
      const SpeciesVCMap& svcmap = *_vcMapVector[m];
      const SpeciesBCMap& sbcmap = *_bcMapVector[m];
      SpeciesFields& sFields = *_speciesFieldsVector[m];
      //MFRPtr& iNorm = *_initialNormVector[m];

      for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        const StorageSite& allFaces = mesh.getFaces();
 
	typename SpeciesVCMap::const_iterator pos = svcmap.find(mesh.getID());
        if (pos==svcmap.end())
	{   
           throw CException("SpeciesModel: Error in VC Map");
	}

	const SpeciesVC<T>& vc = *(pos->second);
     

	//mass fraction
	
        shared_ptr<TArray> mFCell(new TArray(cells.getCount()));
	
	if (!mesh.isDoubleShell())
	  {
	    *mFCell = vc["initialMassFraction"];
	  }
	else
	  {
	    // double shell mesh cells take on values of parents
	    typename SpeciesVCMap::const_iterator posParent = svcmap.find(mesh.getParentMeshID());
	    typename SpeciesVCMap::const_iterator posOther = svcmap.find(mesh.getOtherMeshID());
	    const SpeciesVC<T>& vcP = *(posParent->second);
	    const SpeciesVC<T>& vcO = *(posOther->second);

	    const int NumberOfCells = cells.getCount();
	    for (int c=0; c<NumberOfCells; c++)
	      {
		if ((c < NumberOfCells/4)||(c >= 3*NumberOfCells/4))
		  {
		    (*mFCell)[c] = vcP["initialMassFraction"];
		  }
		else
		  {
		    (*mFCell)[c] = vcO["initialMassFraction"];
		  }
	      }
	  }
	/*
	//Initialize to exp for ln term testing in electric model
	const VectorT3Array& cellCentroid =
	  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

	for (int c=0; c<cells.getCount(); c++)
	{
	  VectorT3 CellCent = cellCentroid[c]; 
	  (*mFCell)[c] = exp(pow(CellCent[0],2));
	  }*/


	/*if (n==0)
	  {
	    *mFCell =_options["initialMassFraction0"];
	  }
	else if (n==1)
	  {
	    *mFCell = _options["initialMassFraction1"];
	  }
	else if (n==2)
	  {
	    *mFCell = _options["initialMassFraction2"];
	  }
	else if (n==3)
	  {
	    *mFCell = _options["initialMassFraction3"];
	  }
	else if (n==4)
	  {
	    *mFCell = _options["initialMassFraction4"];
	  }
	else
	  {
	    *mFCell =_options["initialMassFraction0"];
	    }
	
	if (mesh.isDoubleShell())
	  {
	    char parentInitial[21];
	    char otherInitial[21];
	    sprintf(parentInitial, "initialMassFraction%d", mesh.getParentMeshID());
	    sprintf(otherInitial, "initialMassFraction%d", mesh.getOtherMeshID());

	    for (int c=0; c<cells.getSelfCount(); c++)
	      {
		if (c < cells.getSelfCount()/2)
		  {
		    (*mFCell)[c] = _options[parentInitial];
		  }
		else
		  {
		    (*mFCell)[c] = _options[otherInitial];
		  }
	      }
	  }*/


	//*mFCell = _options["initialMassFraction"];

	/* //Initialize to 1-D linear solution in x direction
	const VectorT3Array& cellCentroid =
	  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

	for (int c=0; c<cells.getCount(); c++)
	{
	  VectorT3 CellCent = cellCentroid[c]; 
	  (*mFCell)[c] = -0.05*CellCent[0] + 0.5;
	  } */

	/*//Initialize to specialized concentrations (for now)
	const VectorT3Array& cellCentroid =
	  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

	for (int c=0; c<cells.getCount(); c++)
	  {
	    VectorT3 CellCent = cellCentroid[c]; 
	    if (n==0)
	      {
		(*mFCell)[c] = 1000.0;
	      }
	    elseif (n==1)
	      {
		(*mFCell)[c] = 0.0;

		(*mFCell)[c] = -0.05*CellCent[0] + 0.5;
		}*/
	
	sFields.massFraction.addArray(cells,mFCell);

	if (_options.ButlerVolmer)
	  {
	    sFields.massFractionElectricModel.addArray(cells,dynamic_pointer_cast<ArrayBase>(mFCell->newCopy()));
	  }

	if (_options.transient)
        {
            sFields.massFractionN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(mFCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
              sFields.massFractionN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(mFCell->newCopy()));

	}

	//diffusivity
        shared_ptr<TArray> diffusivityCell(new TArray(cells.getCount()));
        *diffusivityCell = vc["massDiffusivity"];
        sFields.diffusivity.addArray(cells,diffusivityCell);
	
	//source 
	shared_ptr<TArray> sCell(new TArray(cells.getCount()));
	*sCell = T(0.);
	sFields.source.addArray(cells,sCell);

	//create a zero field
	shared_ptr<TArray> zeroCell(new TArray(cells.getCount()));
	*zeroCell = T(0.0);
	sFields.zero.addArray(cells,zeroCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCount()));
	*oneCell = T(1.0);
	sFields.one.addArray(cells,oneCell);

	//initial temparature gradient array
	shared_ptr<TGradArray> gradT(new TGradArray(cells.getCount()));
	gradT->zero();
	_speciesModelFields.speciesGradient.addArray(cells,gradT);
        
	//inital convection flux at faces
	shared_ptr<TArray> convFlux(new TArray(allFaces.getCount()));
	convFlux->zero();
	sFields.convectionFlux.addArray(allFaces,convFlux);

	//mass flux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            sFields.massFlux.addArray(faces,fluxFace);
          
        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
          
            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            sFields.massFlux.addArray(faces,fluxFace);
          
        }

	//electric potential (for coupling to ElectricModel)
	// only needed once for each mesh (not depenedent on species)
        if (m==0)
	  {
	    shared_ptr<TArray> ePotCell(new TArray(cells.getCount()));
	    ePotCell->zero();
	    sFields.elecPotential.addArray(cells,ePotCell);
	  }
      }	
    
    sFields.diffusivity.syncLocal();
    //iNorm = MFRPtr();
    }
    _niters  =0;    
  }

  SpeciesFields& getSpeciesFields(const int speciesId) {return *_speciesFieldsVector[speciesId];}
  
  SpeciesVCMap& getVCMap(const int speciesId) {return *_vcMapVector[speciesId];}
  SpeciesBCMap& getBCMap(const int speciesId) {return *_bcMapVector[speciesId];}

  SpeciesBC<T>& getBC(const int id, const int speciesId ) {return *(*_bcMapVector[speciesId])[id];}

  SpeciesModelOptions<T>& getOptions() {return _options;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();

    for (int m=0; m<_nSpecies; m++)
    {
        SpeciesFields& sFields = *_speciesFieldsVector[m];

        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];

            const StorageSite& cells = mesh.getCells();
            TArray& mF =
              dynamic_cast<TArray&>(sFields.massFraction[cells]);
            TArray& mFN1 =
              dynamic_cast<TArray&>(sFields.massFractionN1[cells]);

            if (_options.timeDiscretizationOrder > 1)
            {
                TArray& mFN2 =
                  dynamic_cast<TArray&>(sFields.massFractionN2[cells]);
                mFN2 = mFN1;
            }
            mFN1 = mF;
        }
    }
  }

  void initLinearization(LinearSystem& ls, const int& m)
  {
    SpeciesFields& sFields = *_speciesFieldsVector[m];
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&sFields.massFraction,&cells);

        ls.getX().addArray(tIndex,sFields.massFraction.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&sFields.massFlux,&faces);
            ls.getX().addArray(fIndex,sFields.massFlux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&sFields.massFlux,&faces);
            ls.getX().addArray(fIndex,sFields.massFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
     
    }
  }

  void linearize(LinearSystem& ls, const int& m)
  {
    const SpeciesVCMap& svcmap = *_vcMapVector[m];
    const SpeciesBCMap& sbcmap = *_bcMapVector[m];
    SpeciesFields& sFields = *_speciesFieldsVector[m];

    GradientModel<T> speciesGradientModel(_meshes,sFields.massFraction,
					  _speciesModelFields.speciesGradient,_geomFields);
    speciesGradientModel.compute();
    
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  sFields.massFraction,
	  sFields.diffusivity,
	  _speciesModelFields.speciesGradient));
    discretizations.push_back(dd);
    
    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  sFields.massFraction,
	  sFields.convectionFlux,
	  sFields.zero,
	  _speciesModelFields.speciesGradient,
          _options.useCentralDifference));
    discretizations.push_back(cd);
    
    
    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  sFields.massFraction,
	  sFields.source));
    discretizations.push_back(sd);
    
    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativeDiscretization<T,T,T>
             (_meshes,_geomFields,
              sFields.massFraction,
              sFields.massFractionN1,
              sFields.massFractionN2,
              sFields.one,
              _options["timeStep"]));
        
        discretizations.push_back(td);
    }

    shared_ptr<Discretization>
      ibm(new GenericIBDiscretization<T,T,T>
	  (_meshes,_geomFields,sFields.massFraction));
      
    discretizations.push_back(ibm);

    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();

    /* linearize shell mesh */

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	if (mesh.isDoubleShell())
	  {
	    const int parentMeshID = mesh.getParentMeshID();
	    const int otherMeshID = mesh.getOtherMeshID();
	    const Mesh& parentMesh = *_meshes[parentMeshID];
	    const Mesh& otherMesh = *_meshes[otherMeshID];
	  
	    if (_options.ButlerVolmer)
	      {
		bool Cathode = false;
		bool Anode = false;
		if (n == _options["ButlerVolmerCathodeShellMeshID"])
		  {
		    Cathode = true;
		  }
		else if (n == _options["ButlerVolmerAnodeShellMeshID"])
		  {
		    Anode = true;
		  }
		
		LinearizeSpeciesInterface<T, T, T> lbv (_geomFields,
							_options["ButlerVolmerRRConstant"],
							_options["A_coeff"],
							_options["B_coeff"],
							_options["interfaceUnderRelax"],
							Anode,
							Cathode,
							sFields.massFraction,
							sFields.massFractionElectricModel,
							sFields.elecPotential);

		lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );

	      }
	    else
	      {
		LinearizeInterfaceJump<T, T, T> lsm (_options["A_coeff"],
						     _options["B_coeff"],
						     sFields.massFraction);

		lsm.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	      }
	  }
      }

    /* boundary and interface condition */

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            typename SpeciesBCMap::const_iterator pos = sbcmap.find(fg.id);
            if (pos==sbcmap.end())
	    {   
               throw CException("SpeciesModel: Error in BC Map");
	    }

	    const SpeciesBC<T>& bc = *(pos->second);

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  sFields.massFraction,
                                  sFields.massFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedMassFraction")
            {
                FloatValEvaluator<T>
                  bT(bc.getVal("specifiedMassFraction"),faces);
                if (sFields.convectionFlux.hasArray(faces))
                {
                    const TArray& convectingFlux =
                      dynamic_cast<const TArray&>
                      (sFields.convectionFlux[faces]);
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
            else if (bc.bcType == "SpecifiedMassFlux")
            {
                const T specifiedFlux(bc["specifiedMassFlux"]);
                gbc.applyNeumannBC(specifiedFlux);
            }
            else if ((bc.bcType == "Symmetry"))
            {
                 T zeroFlux(NumTypeTraits<T>::getZero());
                 gbc.applyNeumannBC(zeroFlux);
            }
            else
              throw CException(bc.bcType + " not implemented for SpeciesModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  sFields.massFraction,
                                  sFields.massFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }
  
  T getMassFluxIntegral(const Mesh& mesh, const int faceGroupId, const int m)
  {
    SpeciesFields& sFields = *_speciesFieldsVector[m];

    T r(0.);
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const TArray& massFlux =
              dynamic_cast<const TArray&>(sFields.massFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += massFlux[f];
            found=true;
        }
    }

    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
	  {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const TArray& massFlux =
              dynamic_cast<const TArray&>(sFields.massFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += massFlux[f];
            found=true;
	  }
    }


    if (!found)
      throw CException("getMassFluxIntegral: invalid faceGroupID");
    return r;
  }

T getAverageMassFraction(const Mesh& mesh, const int m)
  {
    SpeciesFields& sFields = *_speciesFieldsVector[m];
    const StorageSite& cells = mesh.getCells();
    const TArray& mFCell = dynamic_cast<const TArray&>(sFields.massFraction[cells]);
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    T totalVolume = 0.0;
    T weightedMassFraction = 0.0;
    for(int c=0; c<cells.getSelfCount(); c++)
      {
	totalVolume += cellVolume[c];
	weightedMassFraction += mFCell[c]*cellVolume[c];
      }
    return weightedMassFraction/totalVolume;
  }


  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)
    { 
      bool allConverged=true;
      for (int m=0; m<_nSpecies; m++)
      {
        MFRPtr& iNorm = *_initialNormVector[m];
	MFRPtr& rCurrent = *_currentResidual[m];

        LinearSystem ls;
        initLinearization(ls, m);
        
        ls.initAssembly();

        linearize(ls, m);

        ls.initSolve();

        MFRPtr rNorm(_options.getLinearSolver().solve(ls));

        if (!iNorm) iNorm = rNorm;
        
        MFRPtr normRatio((*rNorm)/(*iNorm));

        cout << "Species Number: " << m << endl;
        cout << _niters << ": " << *rNorm << endl;
        
        _options.getLinearSolver().cleanup();

        ls.postSolve();
        ls.updateSolution();

        _niters++;

	rCurrent = rNorm;
 
        if (!(*rNorm < _options.absoluteTolerance ||
	    *normRatio < _options.relativeTolerance))
	    allConverged=false;
      }
      if (allConverged)
      	break;
    }
  }
    
  void printBCs()
  {
    for (int m=0; m<_nSpecies; m++)
    {
      const SpeciesBCMap& sbcmap = *_bcMapVector[m];
      cout << "Species Number :" << m << endl;

      for(typename SpeciesBCMap::const_iterator pos=sbcmap.begin();
          pos!=sbcmap.end();
          ++pos)
      {
          cout << "Face Group " << pos->first << ":" << endl;
          cout << "    bc type " << pos->second->bcType << endl;
          foreach(typename SpeciesBC<T>::value_type& vp, *(pos->second))
          {
              cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
          }
      }
    }
  }

  T getMassFractionResidual(const int speciesId)
  {
    MFRPtr& rCurrent = *_currentResidual[speciesId];
    SpeciesFields& sFields = *_speciesFieldsVector[speciesId];
    const Field& massFractionField = sFields.massFraction;
    ArrayBase& residualArray = (*rCurrent)[massFractionField];
    T *arrayData = (T*)(residualArray.getData());
    const T residual = arrayData[0];
    return residual;
  }


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  
  vector<SpeciesFields*> _speciesFieldsVector;

  vector<SpeciesBCMap*> _bcMapVector;
  vector<SpeciesVCMap*> _vcMapVector;

  SpeciesModelOptions<T> _options;

  //MFRPtr _initialNorm;
  vector<MFRPtr*> _initialNormVector;
  int _niters;
  const int _nSpecies;

  SpeciesModelFields _speciesModelFields;

  //MFRPtr _currentResidual;
  vector<MFRPtr*> _currentResidual;
};

template<class T>
SpeciesModel<T>::SpeciesModel(const GeomFields& geomFields,
                              const MeshList& meshes,
                              const int nSpecies) :
  Model(meshes),
    _impl(new Impl(geomFields,meshes,nSpecies))
{
  logCtor();
}


template<class T>
SpeciesModel<T>::~SpeciesModel()
{
  logDtor();
}

template<class T>
void
SpeciesModel<T>::init()
{
  _impl->init();
}
  
template<class T>
SpeciesFields&
SpeciesModel<T>::getSpeciesFields(const int speciesId) {return _impl->getSpeciesFields(speciesId);}

template<class T>
typename SpeciesModel<T>::SpeciesBCMap&
SpeciesModel<T>::getBCMap(const int speciesId) {return _impl->getBCMap(speciesId);}

template<class T>
typename SpeciesModel<T>::SpeciesVCMap&
SpeciesModel<T>::getVCMap(const int speciesId) {return _impl->getVCMap(speciesId);}

template<class T>
SpeciesBC<T>&
SpeciesModel<T>::getBC(const int id, const int speciesId) {return _impl->getBC(id, speciesId);}

template<class T>
SpeciesModelOptions<T>&
SpeciesModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
SpeciesModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
void
SpeciesModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}

template<class T>
void
SpeciesModel<T>::updateTime()
{
  _impl->updateTime();
}

template<class T>
T
SpeciesModel<T>::getMassFluxIntegral(const Mesh& mesh, const int faceGroupId, const int m)
{
  return _impl->getMassFluxIntegral(mesh, faceGroupId, m);
}

template<class T>
T
SpeciesModel<T>::getAverageMassFraction(const Mesh& mesh, const int m)
{
  return _impl->getAverageMassFraction(mesh, m);
}

template<class T>
T
SpeciesModel<T>::getMassFractionResidual(const int speciesId)
{
  return _impl->getMassFractionResidual(speciesId);
}
