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
#include "LinearizePotentialInterface.h"
#include "LinearizePCInterface_BV.h"
#include "BatteryElectricDiffusionDiscretization.h"
#include "BatteryPCDiffusionDiscretization.h"
#include "BatteryPC_BCS.h"
#include "BatteryTimeDerivativeDiscretization.h"

template<class T>
class BatteryModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;

  typedef Vector<T,2> VectorT2;
  typedef Array<VectorT2> VectorT2Array;
  typedef SquareTensor<T,2> SquareTensorT2;
  typedef Array<Gradient<VectorT2> > VectorT2GradArray;

  Impl(const GeomFields& geomFields,
       const MeshList& meshes,
       const int nSpecies) :
    _meshes(meshes),
    _geomFields(geomFields),
    _niters(0),
    _nSpecies(nSpecies),
    _batteryModelFields("batteryModel")
  {
    
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

	BatteryPotentialVC<T> *pvc(new BatteryPotentialVC<T>());
        pvc->vcType = "dielectric";
        _pvcMap[mesh.getID()] = pvc;
        
 	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            BatteryPotentialBC<T> *pbc(new BatteryPotentialBC<T>());
            
            _pbcMap[fg.id] = pbc;
            if (fg.groupType == "wall") 
	      {
                pbc->bcType = "SpecifiedPotential";
	      }
	    else if (fg.groupType == "symmetry") 
	      {
                pbc->bcType = "Symmetry";
	      }
	    else
              throw CException("ElectricModel: unknown face group type "
                               + fg.groupType);
	  }
      }

    for (int m=0; m<_nSpecies; m++)
    {
      BatterySpeciesVCMap *svcmap = new BatterySpeciesVCMap();
      _svcMapVector.push_back(svcmap);
      BatterySpeciesBCMap *sbcmap = new BatterySpeciesBCMap();
      _sbcMapVector.push_back(sbcmap);
      BatterySpeciesFields *sFields = new BatterySpeciesFields("species");
      _speciesFieldsVector.push_back(sFields);
      MFRPtr *iNorm = new MFRPtr();
      _initialSpeciesNormVector.push_back(iNorm);
      MFRPtr *rCurrent = new MFRPtr();
      _currentSpeciesResidual.push_back(rCurrent); 
        
      for (int n=0; n<numMeshes; n++)
      {

        const Mesh& mesh = *_meshes[n];
        BatterySpeciesVC<T> *svc(new  BatterySpeciesVC<T>());
        svc->vcType = "flow";
        (*svcmap)[mesh.getID()] = svc;
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
             BatterySpeciesBC<T> *sbc(new  BatterySpeciesBC<T>());
            
            (*sbcmap)[fg.id] = sbc;

            if ((fg.groupType == "wall") ||
                (fg.groupType == "symmetry"))
            {
                sbc->bcType = "SpecifiedMassFlux";
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

    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];	
	const BatteryPotentialVC<T>& pvc = *_pvcMap[mesh.getID()];
	const StorageSite& cells = mesh.getCells();
	const StorageSite& faces = mesh.getFaces();
	
	// initialize fields for electrostatics

	const int nCells = cells.getCount();
	
	//initial potential setup
	shared_ptr<TArray> pCell(new TArray(nCells));

	if (!mesh.isDoubleShell())
	{
	    *pCell = pvc["initialPotential"];
	    }
	else
	  {
	    // double shell mesh cells take on values of parents
	    typename BatteryPotentialVCMap::const_iterator posParent = _pvcMap.find(mesh.getParentMeshID());
	    typename BatteryPotentialVCMap::const_iterator posOther = _pvcMap.find(mesh.getOtherMeshID());
	    const BatteryPotentialVC<T>& pvcP = *(posParent->second);
	    const BatteryPotentialVC<T>& pvcO = *(posOther->second);

	    const int NumberOfCells = cells.getCount();
	    for (int c=0; c<NumberOfCells; c++)
	      {
		if ((c < NumberOfCells/4)||(c >= 3*NumberOfCells/4))
		  {
		    (*pCell)[c] = pvcP["initialPotential"];
		  }
		else
		  {
		    (*pCell)[c] = pvcO["initialPotential"];
		  }
	      }
	  }

	_batteryModelFields.potential.addArray(cells,pCell);

	//conductivity setup
	shared_ptr<TArray> condCell(new TArray(nCells));
	*condCell = pvc["conductivity"];
	_batteryModelFields.conductivity.addArray(cells,condCell);
	  
	// species gradient setup
	shared_ptr<TGradArray> gradS(new TGradArray(cells.getCount()));
	gradS->zero();
	_batteryModelFields.speciesGradient.addArray(cells,gradS);

	//potential gradient setup
	shared_ptr<TGradArray> gradp(new TGradArray(nCells));
	gradp->zero();	
	_batteryModelFields.potential_gradient.addArray(cells,gradp);

	//initial potential flux; Note: potential_flux only stored on boundary faces
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	      
	    shared_ptr<TArray> pFlux(new TArray(faces.getCount()));
	    pFlux->zero();
	    _batteryModelFields.potential_flux.addArray(faces,pFlux);
	      
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	      
	    shared_ptr<TArray> pFlux(new TArray(faces.getCount()));
	    pFlux->zero();
	    _batteryModelFields.potential_flux.addArray(faces,pFlux);
	      
	  }

	shared_ptr<TArray> lnLCCell(new TArray(nCells));
	lnLCCell->zero();
	_batteryModelFields.lnLithiumConcentration.addArray(cells,lnLCCell);

	shared_ptr<TGradArray> gradLnLC(new TGradArray(nCells));
	gradLnLC->zero();
	_batteryModelFields.lnLithiumConcentrationGradient.addArray(cells,gradLnLC);

	//set up combined fields for point-coupled solve
	shared_ptr<VectorT2Array> psCell(new VectorT2Array(nCells));
	psCell->zero();
	_batteryModelFields.potentialAndSpecies.addArray(cells,psCell);

	if (_options.transient)
        {
            _batteryModelFields.potentialAndSpeciesN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(psCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
	      {
		throw CException("BatteryModel: Time Discretization Order > 1 not implimented.");
              /*_batteryModelFields.potentialAndSpeciesN2.addArray(cells,
		dynamic_pointer_cast<ArrayBase>(psCell->newCopy()));*/
	      }

	}

	//combined diffusivity field
	shared_ptr<VectorT2Array> psDiffCell(new VectorT2Array(nCells));
	psDiffCell->zero();
	_batteryModelFields.potentialAndSpeciesDiffusivity.addArray(cells,psDiffCell);	

	//combined gradient setup
	shared_ptr<VectorT2GradArray> gradps(new VectorT2GradArray(nCells));
	gradps->zero();	
	_batteryModelFields.potentialAndSpeciesGradient.addArray(cells,gradps);

	//initial combined flux flied; Note: flux only stored on boundary faces
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	      
	    shared_ptr<VectorT2Array> psFlux(new VectorT2Array(faces.getCount()));
	    psFlux->zero();
	    _batteryModelFields.potentialAndSpeciesFlux.addArray(faces,psFlux);
	      
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	      
	    shared_ptr<VectorT2Array> psFlux(new VectorT2Array(faces.getCount()));
	    psFlux->zero();
	    _batteryModelFields.potentialAndSpeciesFlux.addArray(faces,psFlux);
	      
	  }
      }

    for (int m=0; m<_nSpecies; m++)
    {
      const BatterySpeciesVCMap& svcmap = *_svcMapVector[m];
      const BatterySpeciesBCMap& sbcmap = *_sbcMapVector[m];
      BatterySpeciesFields& sFields = *_speciesFieldsVector[m];
      //MFRPtr& iNorm = *_initialNormVector[m];

      for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        const StorageSite& allFaces = mesh.getFaces();
 
	typename BatterySpeciesVCMap::const_iterator pos = svcmap.find(mesh.getID());
        if (pos==svcmap.end())
	{   
           throw CException("BatteryModel: Error in Species VC Map");
	}

	const BatterySpeciesVC<T>& svc = *(pos->second);
     

	//mass fraction
	
        shared_ptr<TArray> mFCell(new TArray(cells.getCount()));
	
	if (!mesh.isDoubleShell())
	  {
	    *mFCell = svc["initialMassFraction"];
	  }
	else
	  {
	    // double shell mesh cells take on values of parents
	    typename BatterySpeciesVCMap::const_iterator posParent = svcmap.find(mesh.getParentMeshID());
	    typename BatterySpeciesVCMap::const_iterator posOther = svcmap.find(mesh.getOtherMeshID());
	    const BatterySpeciesVC<T>& svcP = *(posParent->second);
	    const BatterySpeciesVC<T>& svcO = *(posOther->second);

	    const int NumberOfCells = cells.getCount();
	    for (int c=0; c<NumberOfCells; c++)
	      {
		if ((c < NumberOfCells/4)||(c >= 3*NumberOfCells/4))
		  {
		    (*mFCell)[c] = svcP["initialMassFraction"];
		  }
		else
		  {
		    (*mFCell)[c] = svcO["initialMassFraction"];
		  }
	      }
	  }
	
	sFields.massFraction.addArray(cells,mFCell);

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
        *diffusivityCell = svc["massDiffusivity"];
        sFields.diffusivity.addArray(cells,diffusivityCell);

	//create a one field
	shared_ptr<TArray> oneCell(new TArray(cells.getCount()));
	*oneCell = T(1.0);
	sFields.one.addArray(cells,oneCell);
        
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

 
      }	
    
    sFields.diffusivity.syncLocal();
    }
    _niters  =0;
    _initialPotentialNorm = MFRPtr();
    _currentPotentialResidual = MFRPtr();
    _initialPCNorm = MFRPtr();
    _currentPCResidual = MFRPtr();
    _batteryModelFields.conductivity.syncLocal();
    copyPCDiffusivity();
    _batteryModelFields.potentialAndSpeciesDiffusivity.syncLocal();
  }

  BatterySpeciesFields& getBatterySpeciesFields(const int speciesId) {return *_speciesFieldsVector[speciesId];}
  BatteryModelFields& getBatteryModelFields() {return _batteryModelFields;}
  BatterySpeciesVCMap& getSpeciesVCMap(const int speciesId) {return *_svcMapVector[speciesId];}
  BatterySpeciesBCMap& getSpeciesBCMap(const int speciesId) {return *_sbcMapVector[speciesId];}
  BatteryPotentialVCMap& getPotentialVCMap() {return _pvcMap;}
  BatteryPotentialBCMap& getPotentialBCMap() {return _pbcMap;}

  BatteryModelOptions<T>& getOptions() {return _options;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();

    for (int m=0; m<_nSpecies; m++)
    {
        BatterySpeciesFields& sFields = *_speciesFieldsVector[m];

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

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];

	const StorageSite& cells = mesh.getCells();
	VectorT2Array& pAndS =
	  dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpecies[cells]);
	VectorT2Array& pAndSN1 =
	  dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesN1[cells]);

	if (_options.timeDiscretizationOrder > 1)
	  {
	    VectorT2Array& pAndSN2 =
	      dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesN2[cells]);
	    pAndSN2 = pAndSN1;
	  }
	pAndSN1 = pAndS;
      }
  }
  
  void initSpeciesLinearization(LinearSystem& ls, const int& SpeciesNumber)
  {
    BatterySpeciesFields& sFields = *_speciesFieldsVector[SpeciesNumber];
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

  void initPCLinearization(LinearSystem& ls)
  {
    BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //only point coupling lithium for now
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_batteryModelFields.potentialAndSpecies,&cells);

        ls.getX().addArray(tIndex,_batteryModelFields.potentialAndSpecies.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<SquareTensorT2,SquareTensorT2,VectorT2>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialAndSpeciesFlux,&faces);
            ls.getX().addArray(fIndex,_batteryModelFields.potentialAndSpeciesFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<SquareTensorT2,VectorT2>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<SquareTensorT2,VectorT2>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialAndSpeciesFlux,&faces);
            ls.getX().addArray(fIndex,_batteryModelFields.potentialAndSpeciesFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<SquareTensorT2,VectorT2>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<SquareTensorT2,VectorT2>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
     
    } 
  } 

 void initPotentialLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_batteryModelFields.potential,&cells);

        ls.getX().addArray(tIndex,_batteryModelFields.potential.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_batteryModelFields.potential_flux,&faces);
            ls.getX().addArray(fIndex,_batteryModelFields.potential_flux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&_batteryModelFields.potential_flux,&faces);
            ls.getX().addArray(fIndex,_batteryModelFields.potential_flux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

    }
  }
  
  void linearizeSpecies(LinearSystem& ls, const int& m)
  {
    const BatterySpeciesVCMap& svcmap = *_svcMapVector[m];
    const BatterySpeciesBCMap& sbcmap = *_sbcMapVector[m];
    BatterySpeciesFields& sFields = *_speciesFieldsVector[m];

    GradientModel<T> speciesGradientModel(_meshes,sFields.massFraction,
					  _batteryModelFields.speciesGradient,_geomFields);
    speciesGradientModel.compute();
    
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  sFields.massFraction,
	  sFields.diffusivity,
	  _batteryModelFields.speciesGradient));
    discretizations.push_back(dd);
    
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

    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();

    // linearize shell mesh
    
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
							T(1.0),
							T(0.0),
							_options["interfaceSpeciesUnderRelax"],
							Anode,
							Cathode,
							sFields.massFraction,
							sFields.massFraction,
							_batteryModelFields.potential);

		lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );

	      }
	    else
	      {
		LinearizeInterfaceJump<T, T, T> lsm (T(1.0),
						     T(0.0),
						     sFields.massFraction);

		lsm.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	      }
	  }
      }

    ///////// boundary and interface condition

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            typename BatterySpeciesBCMap::const_iterator pos = sbcmap.find(fg.id);
            if (pos==sbcmap.end())
	    {   
               throw CException("BatteryModel: Error in Species BC Map");
	    }

	    const BatterySpeciesBC<T>& bc = *(pos->second);

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
              throw CException(bc.bcType + " not implemented for Species in BatteryModel");
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
  
void linearizePotential(LinearSystem& ls)
  {

    GradientModel<T> potentialGradientModel(_meshes,_batteryModelFields.potential,
                              _batteryModelFields.potential_gradient,_geomFields);
    potentialGradientModel.compute();
    
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>(_meshes,_geomFields,
                                            _batteryModelFields.potential,
                                            _batteryModelFields.conductivity,
                                            _batteryModelFields.potential_gradient));
    discretizations.push_back(dd);    

    if(_options.ButlerVolmer)
      {
	//populate lnSpeciesConc Field
	for (int n=0; n<_meshes.size(); n++)
	  {
	    const Mesh& mesh = *_meshes[n];
	    const StorageSite& cells = mesh.getCells();
	    BatterySpeciesFields& sFields = *_speciesFieldsVector[0];
	    const TArray& speciesConc = dynamic_cast<const TArray&>(sFields.massFraction[cells]);
	    TArray& lnSpeciesConc = dynamic_cast<TArray&>(_batteryModelFields.lnLithiumConcentration[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		T CellConc = speciesConc[c];
		if (CellConc <= 0.0)
		  {
		    cout << "Error: Cell Concentration <= 0   MeshID: " << n << " CellNum: " << c << endl;
		    CellConc = 0.01;
		  }
		lnSpeciesConc[c] = log(CellConc); 
	      }
	  }

	//compute gradient for ln term discretization
	GradientModel<T> lnSpeciesGradientModel(_meshes,_batteryModelFields.lnLithiumConcentration,
						_batteryModelFields.lnLithiumConcentrationGradient,_geomFields);
	lnSpeciesGradientModel.compute();
	
	//discretize ln term (only affects electrolyte mesh)
	shared_ptr<Discretization>
	  bedd(new BatteryElectricDiffusionDiscretization<T,T,T>(_meshes,_geomFields,
						_batteryModelFields.potential,
						_batteryModelFields.conductivity,
						_batteryModelFields.lnLithiumConcentration,
						_batteryModelFields.lnLithiumConcentrationGradient,
						_options["BatteryElectrolyteMeshID"]));
	discretizations.push_back(bedd);    

      }
    
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
    
    

    
    const int numMeshes = _meshes.size();

    /* linearize double shell mesh for interface potential jump*/
    
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

	      BatterySpeciesFields& sFields = *_speciesFieldsVector[0];
	      LinearizePotentialInterface<T, T, T> lbv (_geomFields,
							_batteryModelFields.potential,
							sFields.massFraction,
							_options["ButlerVolmerRRConstant"],
							T(1.0),
							T(0.0),
							Anode,
							Cathode);

	      lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	    }
	  else
	    {
	  
	      LinearizeInterfaceJump<T, T, T> lsm (T(1.0),
						   T(0.0),
						   _batteryModelFields.potential);

	      lsm.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	    }
	}
	}

    /* boundary and interface condition */

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
	
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
            const BatteryPotentialBC<T>& bc = *_pbcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.potential,
                                  _batteryModelFields.potential_flux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedPotential")
            {
            	FloatValEvaluator<T> bT(bc.getVal("specifiedPotential"), faces);
                for(int f=0; f<nFaces; f++)
		{
                    gbc.applyDirichletBC(f, bT[f]);
		}
            }
            else if (bc.bcType == "SpecifiedPotentialFlux")
            {
                const T specifiedFlux(bc["specifiedPotentialFlux"]);
                gbc.applyNeumannBC(specifiedFlux);
            }
	    else if (bc.bcType == "Symmetry")
            {
                T zeroFlux(NumTypeTraits<T>::getZero());
                gbc.applyNeumannBC(zeroFlux);
            }
	    else
                throw CException(bc.bcType + " not implemented for Potential in BatteryModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
	  
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

	    
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.potential,
                                  _batteryModelFields.potential_flux, 
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }

void linearizePC(LinearSystem& ls)
  {
    // only lithium (species 0) for now
    const BatterySpeciesVCMap& svcmap = *_svcMapVector[0];
    const BatterySpeciesBCMap& sbcmap = *_sbcMapVector[0];

    GradientModel<VectorT2> psGradientModel(_meshes,_batteryModelFields.potentialAndSpecies,
					  _batteryModelFields.potentialAndSpeciesGradient,_geomFields);
    psGradientModel.compute();
    
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new BatteryPCDiffusionDiscretization<VectorT2,SquareTensorT2,SquareTensorT2>
	 (_meshes,_geomFields,
	  _batteryModelFields.potentialAndSpecies,
	  _batteryModelFields.potentialAndSpeciesDiffusivity,
	  _batteryModelFields.potentialAndSpeciesGradient));
    discretizations.push_back(dd);
    
    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new BatteryTimeDerivativeDiscretization<VectorT2,SquareTensorT2,SquareTensorT2>
             (_meshes,_geomFields,
              _batteryModelFields.potentialAndSpecies,
              _batteryModelFields.potentialAndSpeciesN1,
              _batteryModelFields.potentialAndSpeciesN2,
              _options["timeStep"]));
        
        discretizations.push_back(td);
	}

    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();

    // linearize shell mesh
    
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
		
		LinearizePCInterface_BV<VectorT2, SquareTensorT2, SquareTensorT2> lbv (_geomFields,
							_options["ButlerVolmerRRConstant"],
							_options["interfaceSpeciesUnderRelax"],
							Anode,
							Cathode,
							_batteryModelFields.potentialAndSpecies);

							lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );

	      }
	    else
	      {
		LinearizeInterfaceJump<VectorT2, SquareTensorT2, SquareTensorT2> lsm (T(1.0),
										      NumTypeTraits<VectorT2>::getZero(),
						     _batteryModelFields.potentialAndSpecies);

		lsm.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	      }
	  }
      }

    ///////// boundary and interface condition
    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            typename BatterySpeciesBCMap::const_iterator pos = sbcmap.find(fg.id);
            if (pos==sbcmap.end())
	    {   
               throw CException("BatteryModel: Error in Species BC Map");
	    }

	    const BatterySpeciesBC<T>& sbc = *(pos->second);

	    const BatteryPotentialBC<T>& pbc = *_pbcMap[fg.id];

            BatteryPC_BCS<VectorT2,SquareTensorT2,SquareTensorT2> bbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.potentialAndSpecies,
                                  _batteryModelFields.potentialAndSpeciesFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (sbc.bcType == "SpecifiedMassFraction")
            {
                FloatValEvaluator<T>
                  sbT(sbc.getVal("specifiedMassFraction"),faces);
    
		bbc.applySingleEquationDirichletBC(sbT,1);
            }
            else if (sbc.bcType == "SpecifiedMassFlux")
            {
	      const T specifiedFlux(sbc["specifiedMassFlux"]);
	      bbc.applySingleEquationNeumannBC(specifiedFlux,1);
            }
            else if ((sbc.bcType == "Symmetry"))
            {
	      T zeroFlux(NumTypeTraits<T>::getZero());
	      bbc.applySingleEquationNeumannBC(zeroFlux,1);
            }
            else
             throw CException(sbc.bcType + " not implemented for Species in BatteryModel");

	    if (pbc.bcType == "SpecifiedPotential")
            {
                FloatValEvaluator<T> pbT(pbc.getVal("specifiedPotential"),faces);
    
		bbc.applySingleEquationDirichletBC(pbT,0);
            }
            else if (pbc.bcType == "SpecifiedPotentialFlux")
            {
	      const T specifiedFlux(pbc["specifiedPotentialFlux"]);
	      bbc.applySingleEquationNeumannBC(specifiedFlux,0);
            }
            else if ((pbc.bcType == "Symmetry"))
            {
	      T zeroFlux(NumTypeTraits<T>::getZero());
	      bbc.applySingleEquationNeumannBC(zeroFlux,0);
            }
            else
             throw CException(pbc.bcType + " not implemented for potential in BatteryModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<VectorT2,SquareTensorT2,SquareTensorT2> gbc(faces,mesh,
                                  _geomFields,
				  _batteryModelFields.potentialAndSpecies,
                                  _batteryModelFields.potentialAndSpeciesFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
	}
    }
  }
  
  
  T getMassFluxIntegral(const Mesh& mesh, const int faceGroupId, const int m)
  {
    BatterySpeciesFields& sFields = *_speciesFieldsVector[m];

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

T getPotentialFluxIntegral(const Mesh& mesh, const int faceGroupId)
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
            const TArray& potentialFlux =
              dynamic_cast<const TArray&>(_batteryModelFields.potential_flux[faces]);
            for(int f=0; f<nFaces; f++)
              r += potentialFlux[f];
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
            const TArray& potentialFlux =
              dynamic_cast<const TArray&>(_batteryModelFields.potential_flux[faces]);
            for(int f=0; f<nFaces; f++)
              r += potentialFlux[f];
            found=true;
	  }
    }


    if (!found)
      throw CException("getPotentialFluxIntegral: invalid faceGroupID");
    return r;
  }

T getAverageMassFraction(const Mesh& mesh, const int m)
  {
    BatterySpeciesFields& sFields = *_speciesFieldsVector[m];
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


  void advanceSpecies(const int niter)
  { 
    for(int n=0; n<niter; n++)
    { 
      bool allConverged=true;
      for (int m=0; m<_nSpecies; m++)
      {
        MFRPtr& iNorm = *_initialSpeciesNormVector[m];
	MFRPtr& rCurrent = *_currentSpeciesResidual[m];

        LinearSystem ls;
        initSpeciesLinearization(ls, m);
        
        ls.initAssembly();

        linearizeSpecies(ls, m);

        ls.initSolve();

        MFRPtr rNorm(_options.getLinearSolverSpecies().solve(ls));

        if (!iNorm) iNorm = rNorm;        
        MFRPtr normRatio((*rNorm)/(*iNorm));

        cout << "Species Number: " << m << endl;
        cout << _niters << ": " << *rNorm << endl;
        
        _options.getLinearSolverSpecies().cleanup();

        ls.postSolve();
        ls.updateSolution();

        _niters++;

	rCurrent = rNorm;
 
        if (!(*rNorm < _options.absoluteSpeciesTolerance ||
	    *normRatio < _options.relativeSpeciesTolerance))
	    allConverged=false;
      }
      if (allConverged)
      	break;
	} 
  }

void advancePotential(const int niter)
  {
    for(int n=0; n<niter; n++)
    { 
        LinearSystem ls;

        initPotentialLinearization(ls);
        
        ls.initAssembly();

        linearizePotential(ls);

        ls.initSolve();

        MFRPtr rNorm(_options.getLinearSolverPotential().solve(ls));

        if (!_initialPotentialNorm) _initialPotentialNorm = rNorm;
        
	_currentPotentialResidual = rNorm;

        MFRPtr normRatio((*rNorm)/(*_initialPotentialNorm));

        cout << "Potential: " << _niters << ": " << *rNorm << endl;
        
        _options.getLinearSolverPotential().cleanup();

        ls.postSolve();
        ls.updateSolution();

        _niters++;

        if (*rNorm < _options.absolutePotentialTolerance ||
            *normRatio < _options.relativePotentialTolerance)
          break;
    }
  }

 void advanceCoupled(const int niter)
 { 
   for(int n=0; n<niter; n++)
     {

       LinearSystem ls;

       initPCLinearization(ls);
        
       ls.initAssembly();

       linearizePC(ls);

       ls.initSolve();

       MFRPtr rNorm = _options.getLinearSolverPC().solve(ls);

       if (!_initialPCNorm) _initialPCNorm = rNorm;

       _currentPCResidual = rNorm;

       MFRPtr normRatio((*rNorm)/(*_initialPCNorm));

       _options.getLinearSolverPC().cleanup();

       ls.postSolve();
    
       ls.updateSolution();

       cout << "Point-Coupled: " << _niters << ": " << *rNorm << endl;

       _niters++;

       if (*rNorm < _options.absolutePCTolerance ||
	   *normRatio < _options.relativePCTolerance)
	 break;
     }
 }

  T getMassFractionResidual(const int speciesId)
  {
    MFRPtr& rCurrent = *_currentSpeciesResidual[speciesId];
    BatterySpeciesFields& sFields = *_speciesFieldsVector[speciesId];
    const Field& massFractionField = sFields.massFraction;
    ArrayBase& residualArray = (*rCurrent)[massFractionField];
    T *arrayData = (T*)(residualArray.getData());
    const T residual = arrayData[0];
    return residual;
  }

  T getPotentialResidual()
  {
    const Field& potentialField = _batteryModelFields.potential;
    ArrayBase& residualArray = (*_currentPotentialResidual)[potentialField];
    T *arrayData = (T*)(residualArray.getData());
    const T residual = arrayData[0];
    return residual;
  }

  T getPCResidual(const int v)
  {
    const Field& potentialAndSpeciesField = _batteryModelFields.potentialAndSpecies;
    ArrayBase& residualArray = (*_currentPCResidual)[potentialAndSpeciesField];
    T *arrayData = (T*)(residualArray.getData());
    const T residual = arrayData[v];
    return residual;
  }

  void copySeparateToCoupled()
  {
    for (int n=0; n<_meshes.size(); n++)
      {
	//copy massFraction and potential
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium
	const TArray& mFSeparate = dynamic_cast<const TArray&>(sFields.massFraction[cells]);
	const TArray& pSeparate = dynamic_cast<const TArray&>(_batteryModelFields.potential[cells]);
	VectorT2Array& coupledValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpecies[cells]);

	for (int c=0; c<cells.getCount(); c++)
	  {
	    T massFraction = mFSeparate[c];
	    T potential = pSeparate[c];
	    VectorT2& coupledCellVector = coupledValues[c];

	    //populate coupled field with separate field data
	    coupledCellVector[0] = potential;
	    coupledCellVector[1] = massFraction;
	  }

	//copy fluxes
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const TArray& mFFluxSeparate = dynamic_cast<const TArray&>(sFields.massFlux[faces]);
	    const TArray& pFluxSeparate = dynamic_cast<const TArray&>(_batteryModelFields.potential_flux[faces]);  
	    VectorT2Array& coupledFluxValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesFlux[faces]);
	    for (int f=0; f<faces.getCount(); f++)
	      {
		T massFractionFlux = mFFluxSeparate[f];
		T potentialFlux = pFluxSeparate[f];
		VectorT2& coupledCellFluxVector = coupledFluxValues[f];

		//populate coupled field with separate field data
		coupledCellFluxVector[0] = potentialFlux;
		coupledCellFluxVector[1] = massFractionFlux;
	      }	      
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const TArray& mFFluxSeparate = dynamic_cast<const TArray&>(sFields.massFlux[faces]);
	    const TArray& pFluxSeparate = dynamic_cast<const TArray&>(_batteryModelFields.potential_flux[faces]);  
	    VectorT2Array& coupledFluxValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesFlux[faces]);
	    for (int f=0; f<faces.getCount(); f++)
	      {
		T massFractionFlux = mFFluxSeparate[f];
		T potentialFlux = pFluxSeparate[f];
		VectorT2& coupledCellFluxVector = coupledFluxValues[f];

		//populate coupled field with separate field data
		coupledCellFluxVector[0] = potentialFlux;
		coupledCellFluxVector[1] = massFractionFlux;
	      }	       
	  }

	// copy to N1 and N2 if transient
	if (_options.transient)
	  {
	    const VectorT2Array& coupledValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialAndSpecies[cells]);
	    VectorT2Array& coupledValuesN1 = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesN1[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		coupledValuesN1[c] = coupledValues[c];	
	      }
	      
	    if (_options.timeDiscretizationOrder > 1)
	      {
		VectorT2Array& coupledValuesN2 = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesN2[cells]);
		for (int c=0; c<cells.getCount(); c++)
		  {
		    coupledValuesN2[c] = coupledValues[c];	
		  }
	      }
	  }
      }
  }

  void copyCoupledToSeparate()
  {
    for (int n=0; n<_meshes.size(); n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium

	// copy mass Fraction and potential
	TArray& mFSeparate = dynamic_cast<TArray&>(sFields.massFraction[cells]);
	TArray& pSeparate = dynamic_cast<TArray&>(_batteryModelFields.potential[cells]);
	const VectorT2Array& coupledValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialAndSpecies[cells]);
	for (int c=0; c<cells.getCount(); c++)
	  {
	    VectorT2 coupledCellVector = coupledValues[c];
	    T potential = coupledCellVector[0];
	    T massFraction = coupledCellVector[1];
	    //populate separate fields with coupled field data
	    pSeparate[c] = potential;
	    mFSeparate[c] = massFraction;
	  }

	// copy fluxes
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    TArray& mFFluxSeparate = dynamic_cast<TArray&>(sFields.massFlux[faces]);
	    TArray& pFluxSeparate = dynamic_cast<TArray&>(_batteryModelFields.potential_flux[faces]);  
	    const VectorT2Array& coupledFluxValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialAndSpeciesFlux[faces]);
	    for (int f=0; f<faces.getCount(); f++)
	      {
		VectorT2 coupledCellFluxVector = coupledFluxValues[f];
		T potentialFlux = coupledCellFluxVector[0];
		T massFractionFlux = coupledCellFluxVector[1];

		//populate separate fields with coupled field data
		pFluxSeparate[f] = potentialFlux;
		mFFluxSeparate[f] = massFractionFlux;
	      }	      
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    TArray& mFFluxSeparate = dynamic_cast<TArray&>(sFields.massFlux[faces]);
	    TArray& pFluxSeparate = dynamic_cast<TArray&>(_batteryModelFields.potential_flux[faces]);  
	    const VectorT2Array& coupledFluxValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialAndSpeciesFlux[faces]);
	    for (int f=0; f<faces.getCount(); f++)
	      {
		VectorT2 coupledCellFluxVector = coupledFluxValues[f];
		T potentialFlux = coupledCellFluxVector[0];
		T massFractionFlux = coupledCellFluxVector[1];

		//populate separate fields with coupled field data
		pFluxSeparate[f] = potentialFlux;
		mFFluxSeparate[f] = massFractionFlux;
	      }	      
	  }
      }
  }

void copyPCDiffusivity()
  {
    for (int n=0; n<_meshes.size(); n++)
      {
	//copy massFraction diff and potential cond
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium
	const TArray& mFDiffSeparate = dynamic_cast<const TArray&>(sFields.diffusivity[cells]);
	const TArray& pDiffSeparate = dynamic_cast<const TArray&>(_batteryModelFields.conductivity[cells]);
	VectorT2Array& coupledDiffValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialAndSpeciesDiffusivity[cells]);

	for (int c=0; c<cells.getCount(); c++)
	  {
	    T massDiff = mFDiffSeparate[c];
	    T potentialDiff = pDiffSeparate[c];
	    VectorT2& coupledCellDiffVector = coupledDiffValues[c];

	    //populate coupled field with separate field data
	    coupledCellDiffVector[0] = potentialDiff;
	    coupledCellDiffVector[1] = massDiff;
	  }
      }
  }

private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  
  vector<BatterySpeciesFields*> _speciesFieldsVector;

  vector<BatterySpeciesBCMap*> _sbcMapVector;
  vector<BatterySpeciesVCMap*> _svcMapVector;
  BatteryPotentialBCMap _pbcMap;
  BatteryPotentialVCMap _pvcMap;

  BatteryModelOptions<T> _options;

  MFRPtr _initialPotentialNorm;
  MFRPtr _initialPCNorm;
  vector<MFRPtr*> _initialSpeciesNormVector;
  int _niters;
  const int _nSpecies;

  BatteryModelFields _batteryModelFields;

  MFRPtr _currentPotentialResidual;
  MFRPtr _currentPCResidual;
  vector<MFRPtr*> _currentSpeciesResidual;
};

template<class T>
BatteryModel<T>::BatteryModel(const GeomFields& geomFields,
                              const MeshList& meshes,
                              const int nSpecies) :
  Model(meshes),
    _impl(new Impl(geomFields,meshes,nSpecies))
{
  logCtor();
}


template<class T>
BatteryModel<T>::~BatteryModel()
{
  logDtor();
}

template<class T>
void
BatteryModel<T>::init()
{
  _impl->init();
}
template<class T>
BatterySpeciesFields&
BatteryModel<T>::getBatterySpeciesFields(const int speciesId) {return _impl->getBatterySpeciesFields(speciesId);}

template<class T>
typename BatteryModel<T>::BatterySpeciesBCMap&
BatteryModel<T>::getSpeciesBCMap(const int speciesId) {return _impl->getSpeciesBCMap(speciesId);}

template<class T>
typename BatteryModel<T>::BatterySpeciesVCMap&
BatteryModel<T>::getSpeciesVCMap(const int speciesId) {return _impl->getSpeciesVCMap(speciesId);}

template<class T>
BatteryModelFields&
BatteryModel<T>::getBatteryModelFields() {return _impl->getBatteryModelFields();}

template<class T>
typename BatteryModel<T>::BatteryPotentialBCMap&
BatteryModel<T>::getPotentialBCMap() {return _impl->getPotentialBCMap();}

template<class T>
typename BatteryModel<T>::BatteryPotentialVCMap&
BatteryModel<T>::getPotentialVCMap() {return _impl->getPotentialVCMap();}

template<class T>
BatteryModelOptions<T>&
BatteryModel<T>::getOptions() {return _impl->getOptions();}

template<class T>
void
BatteryModel<T>::advanceSpecies(const int niter)
{
  _impl->advanceSpecies(niter);
}

template<class T>
void
BatteryModel<T>::advancePotential(const int niter)
{
  _impl->advancePotential(niter);
}

template<class T>
void
BatteryModel<T>::advanceCoupled(const int niter)
{
  _impl->advanceCoupled(niter);
}

template<class T>
void
BatteryModel<T>::updateTime()
{
  _impl->updateTime();
}

template<class T>
T
BatteryModel<T>::getMassFluxIntegral(const Mesh& mesh, const int faceGroupId, const int m)
{
  return _impl->getMassFluxIntegral(mesh, faceGroupId, m);
}

template<class T>
T
BatteryModel<T>::getPotentialFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
  return _impl->getPotentialFluxIntegral(mesh, faceGroupId);
}

template<class T>
T
BatteryModel<T>::getAverageMassFraction(const Mesh& mesh, const int m)
{
  return _impl->getAverageMassFraction(mesh, m);
}

template<class T>
T
BatteryModel<T>::getMassFractionResidual(const int speciesId)
{
  return _impl->getMassFractionResidual(speciesId);
}

template<class T>
T
BatteryModel<T>::getPotentialResidual()
{
  return _impl->getPotentialResidual();
}

template<class T>
T
BatteryModel<T>::getPCResidual(const int v)
{
  return _impl->getPCResidual(v);
}

template<class T>
void
BatteryModel<T>::copyCoupledToSeparate()
{
  _impl->copyCoupledToSeparate();
}

template<class T>
void
BatteryModel<T>::copySeparateToCoupled()
{
  _impl->copySeparateToCoupled();
}
