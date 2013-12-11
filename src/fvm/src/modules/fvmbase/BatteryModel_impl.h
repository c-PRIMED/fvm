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
#include "BatteryLinearizeSpeciesInterface.h"
#include "BatteryLinearizePotentialInterface.h"
#include "BatteryLinearizeThermalInterface.h"
#include "BatteryPCLinearizeInterface_BV.h"
#include "BatteryBinaryElectrolyteDiscretization.h"
#include "BatteryPCDiffusionDiscretization.h"
#include "BatteryPC_BCS.h"
#include "BatteryPCTimeDerivativeDiscretization.h"
#include "BatteryPCBinaryElectrolyteDiscretization.h"
#include "BatteryPCHeatSourceDiscretization.h"
#include "BatteryFixInterfaceGhost.h"
#include "LinearizeInterfaceJumpUnconnected.h"

template<class T>
class BatteryModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Gradient<T> TGradType;
  typedef Array<Gradient<T> > TGradArray;
  typedef CRMatrix<T,T,T> T_Matrix;

  // PC Thermal included
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef SquareTensor<T,3> SquareTensorT3;
  typedef Array<Gradient<VectorT3> > VectorT3GradArray;
  typedef Gradient<VectorT3> VectorT3Grad;
  typedef DiagonalTensor<T,3> DiagTensorT3;

  // PC No Thermal included
  typedef Vector<T,2> VectorT2;
  typedef Array<VectorT2> VectorT2Array;
  typedef SquareTensor<T,2> SquareTensorT2;
  typedef Array<Gradient<VectorT2> > VectorT2GradArray;
  typedef Gradient<VectorT2> VectorT2Grad;
  typedef DiagonalTensor<T,2> DiagTensorT2;

  Impl(GeomFields& geomFields,
       const MeshList& realMeshes,
       const MeshList& meshes,
       const int nSpecies) :
    _realMeshes(realMeshes),
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
                //pbc->bcType = "SpecifiedPotential";
		pbc->bcType = "Symmetry";
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

    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      BatteryThermalVC<T> *tvc(new BatteryThermalVC<T>());
      tvc->vcType = "flow";
      _tvcMap[mesh.getID()] = tvc;
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            BatteryThermalBC<T> *tbc(new BatteryThermalBC<T>());
            
            _tbcMap[fg.id] = tbc;

            if ((fg.groupType == "wall") ||
                (fg.groupType == "symmetry"))
            {
	      //tbc->bcType = "SpecifiedHeatFlux";
	      tbc->bcType = "Symmetry";
            }
            else if ((fg.groupType == "velocity-inlet") ||
                     (fg.groupType == "pressure-outlet"))
            {
                tbc->bcType = "SpecifiedTemperature";
            }
            else
              throw CException("ThermalModel: unknown face group type "
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
	      //sbc->bcType = "SpecifiedMassFlux";
	      sbc->bcType = "Symmetry";
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

    // initialize fields for battery model
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];	
	const BatteryPotentialVC<T>& pvc = *_pvcMap[mesh.getID()];
	const BatteryThermalVC<T>& tvc = *_tvcMap[mesh.getID()];
	const StorageSite& cells = mesh.getCells();

	const int nCells = cells.getCount();
	
	// print cell centroids for test case
	//const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
	VectorT3Array& cellCentroid = dynamic_cast<VectorT3Array&>(_geomFields.coordinate[cells]);

	cout << "Mesh: " << n << endl;
	for (int c=0; c<nCells; c++)
	  {
	    if (((cellCentroid[c])[2] < -3.33)&&((cellCentroid[c])[2] > -3.34))
	      cout << c << ": " << (cellCentroid[c])[0] << " " << (cellCentroid[c])[1] << endl;
	  }

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

	if (_options.transient)
	  {
	    // not actually discretized as an unsteady term
	    // only used to recover values from beginning of given time step
	    // in case a new technique needs to be tried after a solution failure
	    _batteryModelFields.potentialN1.addArray(cells, dynamic_pointer_cast<ArrayBase>(pCell->newCopy()));
	  }

	//initial temperature setup
	shared_ptr<TArray> tCell(new TArray(nCells));

	if (!mesh.isDoubleShell())
	  {
	    *tCell = tvc["initialTemperature"];
	  }
	else
	  {
	    // double shell mesh cells take on values of parents
	    typename BatteryThermalVCMap::const_iterator posParent = _tvcMap.find(mesh.getParentMeshID());
	    typename BatteryThermalVCMap::const_iterator posOther = _tvcMap.find(mesh.getOtherMeshID());
	    const BatteryThermalVC<T>& tvcP = *(posParent->second);
	    const BatteryThermalVC<T>& tvcO = *(posOther->second);

	    const int NumberOfCells = cells.getCount();
	    for (int c=0; c<NumberOfCells; c++)
	      {
		if ((c < NumberOfCells/4)||(c >= 3*NumberOfCells/4))
		  {
		    (*tCell)[c] = tvcP["initialTemperature"];
		  }
		else
		  {
		    (*tCell)[c] = tvcO["initialTemperature"];
		  }
	      }
	  }

	_batteryModelFields.temperature.addArray(cells,tCell);

	if (_options.transient)
	  {
	    _batteryModelFields.temperatureN1.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
	      _batteryModelFields.temperatureN2.addArray(cells, dynamic_pointer_cast<ArrayBase>(tCell->newCopy()));
	  }

	//electric/ionic conductivity setup
	shared_ptr<TArray> condCell(new TArray(nCells));
	*condCell = pvc["conductivity"];
	_batteryModelFields.conductivity.addArray(cells,condCell);

	//thermal conductivity setup
	shared_ptr<TArray> thermCondCell(new TArray(nCells));
	*thermCondCell = tvc["thermalConductivity"];
	_batteryModelFields.thermalConductivity.addArray(cells,thermCondCell);
	  
	// species gradient setup
	shared_ptr<TGradArray> gradS(new TGradArray(cells.getCount()));
	gradS->zero();
	_batteryModelFields.speciesGradient.addArray(cells,gradS);

	//potential gradient setup
	shared_ptr<TGradArray> gradp(new TGradArray(nCells));
	gradp->zero();	
	_batteryModelFields.potential_gradient.addArray(cells,gradp);

	//temperature gradient setup
	shared_ptr<TGradArray> gradt(new TGradArray(nCells));
	gradt->zero();	
	_batteryModelFields.temperatureGradient.addArray(cells,gradt);

	//create rho *specific heat field   rho*Cp
	shared_ptr<TArray> rhoCp(new TArray(nCells));
	*rhoCp = tvc["density"] * tvc["specificHeat"];
	_batteryModelFields.rhoCp.addArray(cells, rhoCp);

	//heat source 
	shared_ptr<TArray> sCell(new TArray(nCells));
	*sCell = T(0.);
	_batteryModelFields.heatSource.addArray(cells,sCell);

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

	//initial heat flux; Note: heatflux only stored on boundary faces
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	      
	    shared_ptr<TArray> tFlux(new TArray(faces.getCount()));
	    tFlux->zero();
	    _batteryModelFields.heatFlux.addArray(faces,tFlux);
	      
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	      
	    shared_ptr<TArray> tFlux(new TArray(faces.getCount()));
	    tFlux->zero();
	    _batteryModelFields.heatFlux.addArray(faces,tFlux);
	      
	  }

	shared_ptr<TArray> lnLCCell(new TArray(nCells));
	lnLCCell->zero();
	_batteryModelFields.lnLithiumConcentration.addArray(cells,lnLCCell);

	shared_ptr<TGradArray> gradLnLC(new TGradArray(nCells));
	gradLnLC->zero();
	_batteryModelFields.lnLithiumConcentrationGradient.addArray(cells,gradLnLC);

	//set up combined fields for point-coupled solve
	if (_options.thermalModelPC)
	  {
	    shared_ptr<VectorT3Array> pstCell(new VectorT3Array(nCells));
	    pstCell->zero();
	    _batteryModelFields.potentialSpeciesTemp.addArray(cells,pstCell);

	    if (_options.transient)
	      {
		_batteryModelFields.potentialSpeciesTempN1.addArray(cells,
								    dynamic_pointer_cast<ArrayBase>(pstCell->newCopy()));
		if (_options.timeDiscretizationOrder > 1)
		  {
		    //throw CException("BatteryModel: Time Discretization Order > 1 not implimented.");
		    _batteryModelFields.potentialSpeciesTempN2.addArray(cells,
		      dynamic_pointer_cast<ArrayBase>(pstCell->newCopy()));
		  }

	      }

	    //combined diffusivity field
	    shared_ptr<VectorT3Array> pstDiffCell(new VectorT3Array(nCells));
	    pstDiffCell->zero();
	    _batteryModelFields.potentialSpeciesTempDiffusivity.addArray(cells,pstDiffCell);	

	    //combined gradient setup
	    shared_ptr<VectorT3GradArray> gradpst(new VectorT3GradArray(nCells));
	    gradpst->zero();	
	    _batteryModelFields.potentialSpeciesTempGradient.addArray(cells,gradpst);

	    //initial combined flux flied; Note: flux only stored on boundary faces
	    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
	      
		shared_ptr<VectorT3Array> pstFlux(new VectorT3Array(faces.getCount()));
		pstFlux->zero();
		_batteryModelFields.potentialSpeciesTempFlux.addArray(faces,pstFlux);
	      
	      }
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
	      
		shared_ptr<VectorT3Array> pstFlux(new VectorT3Array(faces.getCount()));
		pstFlux->zero();
		_batteryModelFields.potentialSpeciesTempFlux.addArray(faces,pstFlux);
	      }
	  }
	else
	  {
	    shared_ptr<VectorT2Array> pstCell(new VectorT2Array(nCells));
	    pstCell->zero();
	    _batteryModelFields.potentialSpeciesTemp.addArray(cells,pstCell);

	    if (_options.transient)
	      {
		_batteryModelFields.potentialSpeciesTempN1.addArray(cells,
								    dynamic_pointer_cast<ArrayBase>(pstCell->newCopy()));
		if (_options.timeDiscretizationOrder > 1)
		  {
		    //throw CException("BatteryModel: Time Discretization Order > 1 not implimented.");
		    _batteryModelFields.potentialSpeciesTempN2.addArray(cells,
		      dynamic_pointer_cast<ArrayBase>(pstCell->newCopy()));
		  }

	      }

	    //combined diffusivity field
	    shared_ptr<VectorT2Array> pstDiffCell(new VectorT2Array(nCells));
	    pstDiffCell->zero();
	    _batteryModelFields.potentialSpeciesTempDiffusivity.addArray(cells,pstDiffCell);	

	    //combined gradient setup
	    shared_ptr<VectorT2GradArray> gradpst(new VectorT2GradArray(nCells));
	    gradpst->zero();	
	    _batteryModelFields.potentialSpeciesTempGradient.addArray(cells,gradpst);

	    //initial combined flux flied; Note: flux only stored on boundary faces
	    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
	      
		shared_ptr<VectorT2Array> pstFlux(new VectorT2Array(faces.getCount()));
		pstFlux->zero();
		_batteryModelFields.potentialSpeciesTempFlux.addArray(faces,pstFlux);
	      
	      }
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
	      
		shared_ptr<VectorT2Array> pstFlux(new VectorT2Array(faces.getCount()));
		pstFlux->zero();
		_batteryModelFields.potentialSpeciesTempFlux.addArray(faces,pstFlux);
	      }
	  }
      }

    for (int m=0; m<_nSpecies; m++)
      {
	const BatterySpeciesVCMap& svcmap = *_svcMapVector[m];
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
     

	//concentration 
	
        shared_ptr<TArray> concCell(new TArray(cells.getCount()));
	
	if (!mesh.isDoubleShell())
	  {
	    *concCell = svc["initialConcentration"];
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
		    (*concCell)[c] = svcP["initialConcentration"];
		  }
		else
		  {
		    (*concCell)[c] = svcO["initialConcentration"];
		  }
	      }
	  }
	
	sFields.concentration.addArray(cells,concCell);

	if (_options.transient)
	  {
            sFields.concentrationN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(concCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
              sFields.concentrationN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(concCell->newCopy()));

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
    _initialThermalNorm = MFRPtr();
    _currentThermalResidual = MFRPtr();
    _initialPCNorm = MFRPtr();
    _currentPCResidual = MFRPtr();
    _batteryModelFields.conductivity.syncLocal();
    _batteryModelFields.thermalConductivity.syncLocal();
    copyPCDiffusivity();
    _batteryModelFields.potentialSpeciesTempDiffusivity.syncLocal();
  }

  BatterySpeciesFields& getBatterySpeciesFields(const int speciesId) {return *_speciesFieldsVector[speciesId];}
  BatteryModelFields& getBatteryModelFields() {return _batteryModelFields;}
  BatterySpeciesVCMap& getSpeciesVCMap(const int speciesId) {return *_svcMapVector[speciesId];}
  BatterySpeciesBCMap& getSpeciesBCMap(const int speciesId) {return *_sbcMapVector[speciesId];}
  BatteryPotentialVCMap& getPotentialVCMap() {return _pvcMap;}
  BatteryPotentialBCMap& getPotentialBCMap() {return _pbcMap;}
  BatteryThermalVCMap& getThermalVCMap() {return _tvcMap;}
  BatteryThermalBCMap& getThermalBCMap() {return _tbcMap;}

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
              dynamic_cast<TArray&>(sFields.concentration[cells]);
            TArray& mFN1 =
              dynamic_cast<TArray&>(sFields.concentrationN1[cells]);

            if (_options.timeDiscretizationOrder > 1)
            {
                TArray& mFN2 =
                  dynamic_cast<TArray&>(sFields.concentrationN2[cells]);
                mFN2 = mFN1;
            }
            mFN1 = mF;
        }
    }

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();

	if (_options.thermalModelPC)
	  {
	    VectorT3Array& pST =
	      dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    VectorT3Array& pSTN1 =
	      dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);

	    if (_options.timeDiscretizationOrder > 1)
	      {
		//throw CException("BatteryModel: Time Discretization Order > 1 not implimented.");
		VectorT3Array& pSTN2 = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		pSTN2 = pSTN1;
	      }
	    pSTN1 = pST;
	  }
	else
	  {
	    VectorT2Array& pST =
	      dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    VectorT2Array& pSTN1 =
	      dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);

	    if (_options.timeDiscretizationOrder > 1)
	      {
		//throw CException("BatteryModel: Time Discretization Order > 1 not implimented.");
		VectorT2Array& pSTN2 = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		pSTN2 = pSTN1;
	      }
	    pSTN1 = pST;
	  }

	TArray& temp =
              dynamic_cast<TArray&>(_batteryModelFields.temperature[cells]);
        TArray& tempN1 =
              dynamic_cast<TArray&>(_batteryModelFields.temperatureN1[cells]);

        if (_options.timeDiscretizationOrder > 1)
	  {
	    TArray& tempN2 =
	      dynamic_cast<TArray&>(_batteryModelFields.temperatureN2[cells]);
	    tempN2 = tempN1;
	  }
	tempN1 = temp;

	// only used for recoverLastTimeStep, no unsteady term in potential equation
	TArray& potential =
              dynamic_cast<TArray&>(_batteryModelFields.potential[cells]);
        TArray& potentialN1 =
              dynamic_cast<TArray&>(_batteryModelFields.potentialN1[cells]);
	potentialN1 = potential;

      }
  }

void recoverLastTimestep()
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
              dynamic_cast<TArray&>(sFields.concentration[cells]);
            TArray& mFN1 =
              dynamic_cast<TArray&>(sFields.concentrationN1[cells]);

            mF = mFN1;

	    /*if (_options.timeDiscretizationOrder > 1)
	      {
		TArray& mFN2 =
                  dynamic_cast<TArray&>(sFields.concentrationN2[cells]);
                mFN1 = mFN2;
		}*/
        }
    }

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();

	if (_options.thermalModelPC)
	  {
	    VectorT3Array& pST =
	      dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    VectorT3Array& pSTN1 =
	      dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);
	    pST = pSTN1;

	    /*if (_options.timeDiscretizationOrder > 1)
	      {
		VectorT3Array& pSTN2 = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		pSTN1 = pSTN2;
		}*/
	  }
	else
	  {
	    VectorT2Array& pST =
	      dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    VectorT2Array& pSTN1 =
	      dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);
	    pST = pSTN1;
	    /*if (_options.timeDiscretizationOrder > 1)
	      {
		VectorT2Array& pSTN2 = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		pSTN1 = pSTN2;
		}*/
	  }

	TArray& temp =
              dynamic_cast<TArray&>(_batteryModelFields.temperature[cells]);
        TArray& tempN1 =
              dynamic_cast<TArray&>(_batteryModelFields.temperatureN1[cells]);
	temp = tempN1;

	/*if (_options.timeDiscretizationOrder > 1)
	  {
	    TArray& tempN2 =
	      dynamic_cast<TArray&>(_batteryModelFields.temperatureN2[cells]);
	    tempN1 = tempN2;
	    }*/

	TArray& potential =
              dynamic_cast<TArray&>(_batteryModelFields.potential[cells]);
        TArray& potentialN1 =
              dynamic_cast<TArray&>(_batteryModelFields.potentialN1[cells]);
	potential = potentialN1;
      }
  }
  
  void initSpeciesLinearization(LinearSystem& ls, LinearSystem& lsShell, const int& SpeciesNumber)
  {
    BatterySpeciesFields& sFields = *_speciesFieldsVector[SpeciesNumber];
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
	
	if ((!(mesh.isDoubleShell()))||((mesh.isDoubleShell())&&(mesh.isConnectedShell())))
	  {
	    const StorageSite& cells = mesh.getCells();
	    MultiField::ArrayIndex tIndex(&sFields.concentration,&cells);

	    ls.getX().addArray(tIndex,sFields.concentration.getArrayPtr(cells));

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
	else
	  {
	    // must be unconnected shell, add to shell LS
	     const StorageSite& cells = mesh.getCells();
	    MultiField::ArrayIndex tIndex(&sFields.concentration,&cells);

	    lsShell.getX().addArray(tIndex,sFields.concentration.getArrayPtr(cells));

	    const CRConnectivity& cellCells = mesh.getCellCells();

	    shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

	    lsShell.getMatrix().addMatrix(tIndex,tIndex,m);
	    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;

		MultiField::ArrayIndex fIndex(&sFields.massFlux,&faces);
		lsShell.getX().addArray(fIndex,sFields.massFlux.getArrayPtr(faces));

		const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
		lsShell.getMatrix().addMatrix(fIndex,tIndex,mft);

		shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
		lsShell.getMatrix().addMatrix(fIndex,fIndex,mff);
	      }

	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;

		MultiField::ArrayIndex fIndex(&sFields.massFlux,&faces);
		lsShell.getX().addArray(fIndex,sFields.massFlux.getArrayPtr(faces));

		const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
		lsShell.getMatrix().addMatrix(fIndex,tIndex,mft);

		shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
		lsShell.getMatrix().addMatrix(fIndex,fIndex,mff);
	      }
	  }
    } 
  } 

  void initPCLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

	if (mesh.isDoubleShell())
	  {

	    const StorageSite& cells = mesh.getCells();
	    MultiField::ArrayIndex tIndex(&_batteryModelFields.potentialSpeciesTemp,&cells);

	    ls.getX().addArray(tIndex,_batteryModelFields.potentialSpeciesTemp.getArrayPtr(cells));

	    const CRConnectivity& cellCells = mesh.getCellCells();

	    if (_options.thermalModelPC)
	      {
		shared_ptr<Matrix> m(new CRMatrix<SquareTensorT3,SquareTensorT3,VectorT3>(cellCells));

		ls.getMatrix().addMatrix(tIndex,tIndex,m);

		foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<SquareTensorT3,VectorT3>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<SquareTensorT3,VectorT3>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }

		foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<SquareTensorT3,VectorT3>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<SquareTensorT3,VectorT3>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }
	      }
	    else
	      {
		shared_ptr<Matrix> m(new CRMatrix<SquareTensorT2,SquareTensorT2,VectorT2>(cellCells));

		ls.getMatrix().addMatrix(tIndex,tIndex,m);

		foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

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

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<SquareTensorT2,VectorT2>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<SquareTensorT2,VectorT2>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }
	      }
	  }
	else
	  {
	    const StorageSite& cells = mesh.getCells();
	    MultiField::ArrayIndex tIndex(&_batteryModelFields.potentialSpeciesTemp,&cells);

	    ls.getX().addArray(tIndex,_batteryModelFields.potentialSpeciesTemp.getArrayPtr(cells));

	    const CRConnectivity& cellCells = mesh.getCellCells();

	    if (_options.thermalModelPC)
	      {
		shared_ptr<Matrix> m(new CRMatrix<DiagTensorT3,DiagTensorT3,VectorT3>(cellCells));

		ls.getMatrix().addMatrix(tIndex,tIndex,m);

		foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<DiagTensorT3,VectorT3>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<DiagTensorT3,VectorT3>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }

		foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<DiagTensorT3,VectorT3>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<DiagTensorT3,VectorT3>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }
	      }
	    else
	      {
		shared_ptr<Matrix> m(new CRMatrix<DiagTensorT2,DiagTensorT2,VectorT2>(cellCells));

		ls.getMatrix().addMatrix(tIndex,tIndex,m);

		foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<DiagTensorT2,VectorT2>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<DiagTensorT2,VectorT2>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }

		foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;

		    MultiField::ArrayIndex fIndex(&_batteryModelFields.potentialSpeciesTempFlux,&faces);
		    ls.getX().addArray(fIndex,_batteryModelFields.potentialSpeciesTempFlux.getArrayPtr(faces));

		    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

		    shared_ptr<Matrix> mft(new FluxJacobianMatrix<DiagTensorT2,VectorT2>(faceCells));
		    ls.getMatrix().addMatrix(fIndex,tIndex,mft);

		    shared_ptr<Matrix> mff(new DiagonalMatrix<DiagTensorT2,VectorT2>(faces.getCount()));
		    ls.getMatrix().addMatrix(fIndex,fIndex,mff);
		  }
	      }
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
  
void initThermalLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_batteryModelFields.temperature,&cells);

        ls.getX().addArray(tIndex,_batteryModelFields.temperature.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_batteryModelFields.heatFlux,&faces);
            ls.getX().addArray(fIndex,_batteryModelFields.heatFlux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&_batteryModelFields.heatFlux,&faces);
            ls.getX().addArray(fIndex,_batteryModelFields.heatFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

    }
  }

  void linearizeSpecies(LinearSystem& ls, LinearSystem& lsShell, const int& m)
  {
    const BatterySpeciesBCMap& sbcmap = *_sbcMapVector[m];
    BatterySpeciesFields& sFields = *_speciesFieldsVector[m];

    sFields.concentration.syncLocal();

    //fix interface if it is an unconnected shell mesh 
    // (tested inside so that multiple interfaces can be handled)
    BatteryFixInterfaceGhost<T, T, T> fig (_meshes,
					   _geomFields,
					   sFields.concentration,
					   sFields.diffusivity);

    fig.fixInterfaces();


    GradientModel<T> speciesGradientModel(_meshes,sFields.concentration,
					  _batteryModelFields.speciesGradient,_geomFields);
    speciesGradientModel.compute();
    
    DiscrList discretizations;

    /*shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  sFields.concentration,
	  sFields.diffusivity,
	  _batteryModelFields.speciesGradient));
    discretizations.push_back(dd);
    
    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativeDiscretization<T,T,T>
             (_meshes,_geomFields,
              sFields.concentration,
              sFields.concentrationN1,
              sFields.concentrationN2,
              sFields.one,
              _options["timeStep"]));
        
        discretizations.push_back(td);
    }

      
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
    */

    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_realMeshes,_geomFields,
	  sFields.concentration,
	  sFields.diffusivity,
	  _batteryModelFields.speciesGradient));
    discretizations.push_back(dd);
    
    if (_options.transient)
      {
        shared_ptr<Discretization>
          td(new TimeDerivativeDiscretization<T,T,T>
             (_realMeshes,_geomFields,
              sFields.concentration,
              sFields.concentrationN1,
              sFields.concentrationN2,
              sFields.one,
              _options["timeStep"]));
        
        discretizations.push_back(td);
      }

      
    Linearizer linearizer;

    linearizer.linearize(discretizations,_realMeshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

    const int numMeshes = _meshes.size();

    // linearize shell mesh
    
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	if ((mesh.isDoubleShell())&&(mesh.isConnectedShell()))
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
		
		BatteryLinearizeSpeciesInterface<T, T, T> lbv (_geomFields,
							_options["ButlerVolmerRRConstant"],
							_options["interfaceSpeciesUnderRelax"],
							Anode,
							Cathode,
							sFields.concentration,
							sFields.concentration,
							_batteryModelFields.potential,
							_batteryModelFields.temperature);

		lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );

	      }
	    else
	      {
		LinearizeInterfaceJump<T, T, T> lsm (T(1.0),
						     T(0.0),
						     sFields.concentration);

		lsm.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	      }
	  }
	/*
	if (n==0) 
	  {
	  cout << "BEFORE" << endl;
	  printMatrixElementsOnFace(mesh, 8, ls, sFields.concentration);
	  }
	else if (n==1)
	  printMatrixElementsOnFace(mesh, 11, ls, sFields.concentration);
	*/
	if ((mesh.isDoubleShell())&&(!(mesh.isConnectedShell())))
	  {
	    const int parentMeshID = mesh.getParentMeshID();
	    const int otherMeshID = mesh.getOtherMeshID();
	    const Mesh& parentMesh = *_meshes[parentMeshID];
	    const Mesh& otherMesh = *_meshes[otherMeshID];
	    LinearizeInterfaceJumpUnconnected<T, T, T> lsm (T(1.0),
						     T(0.0),
						     sFields.concentration);

	    lsm.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() , lsShell);

	    //sync for debugging
	    //sync doesn't effect result, but matches what the matrix/residual now contains
	    sFields.concentration.syncLocal();
	  }
      }
    /*
    cout << "AFTER" << endl;
    printMatrixElementsOnFace(*_meshes[0], 8, ls, sFields.concentration);
    printMatrixElementsOnFace(*_meshes[1], 11, ls, sFields.concentration);
    */

    ///////// boundary and interface condition
    const int numRealMeshes = _realMeshes.size();
    for (int n=0; n<numRealMeshes; n++)
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
                                  sFields.concentration,
                                  sFields.massFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedConcentration")
            {
                FloatValEvaluator<T>
                  bT(bc.getVal("specifiedConcentration"),faces);
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
				  sFields.concentration,
				  sFields.massFlux,
				  ls.getMatrix(), ls.getX(), ls.getB());

	    gbc.applyInterfaceBC();
	  }
	  
    }
  }
  
void linearizePotential(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();

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
	bool foundElectrolyte = false;
	//populate lnSpeciesConc Field
	for (int n=0; n<numMeshes; n++)
	  {
	    if (n==_options["BatteryElectrolyteMeshID"])
	      {
		foundElectrolyte = true;
		const Mesh& mesh = *_meshes[n];
		const StorageSite& cells = mesh.getCells();
		BatterySpeciesFields& sFields = *_speciesFieldsVector[0];
		const TArray& speciesConc = dynamic_cast<const TArray&>(sFields.concentration[cells]);
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
	  }

	if ((!(foundElectrolyte))&&(_options.ButlerVolmer))
	  cout << "Warning: Electrolyte Mesh ID not set." << endl;

	//compute gradient for ln term discretization
	GradientModel<T> lnSpeciesGradientModel(_meshes,_batteryModelFields.lnLithiumConcentration,
						_batteryModelFields.lnLithiumConcentrationGradient,_geomFields);
	lnSpeciesGradientModel.compute();
	
	//discretize ln term (only affects electrolyte mesh)
	shared_ptr<Discretization>
	  bedd(new BatteryBinaryElectrolyteDiscretization<T,T,T>(_meshes,_geomFields,
						_batteryModelFields.potential,
						_batteryModelFields.conductivity,
						_batteryModelFields.lnLithiumConcentration,
						_batteryModelFields.lnLithiumConcentrationGradient,
						_batteryModelFields.temperature,
						_options["BatteryElectrolyteMeshID"]));
	discretizations.push_back(bedd);    

      }
    
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
    
    

    
   

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
	      BatteryLinearizePotentialInterface<T, T, T> lbv (_geomFields,
							_batteryModelFields.potential,
							sFields.concentration,
							_batteryModelFields.temperature,
							_options["ButlerVolmerRRConstant"],
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

void linearizeThermal(LinearSystem& ls)
  {

    const int numMeshes = _meshes.size();

    // calculate new thermal gradient
    GradientModel<T> thermalGradientModel(_meshes,_batteryModelFields.temperature,
                              _batteryModelFields.temperatureGradient,_geomFields);
    thermalGradientModel.compute();
    
    // recalculate other gradients for use in heat source term
    GradientModel<T> potentialGradientModel(_meshes,_batteryModelFields.potential,
                              _batteryModelFields.potential_gradient,_geomFields);
    potentialGradientModel.compute();

    BatterySpeciesFields& sFields = *_speciesFieldsVector[0];
    GradientModel<T> speciesGradientModel(_meshes,sFields.concentration,
					  _batteryModelFields.speciesGradient,_geomFields);
    speciesGradientModel.compute();
    
   
    // fill source field with dot product of gradients from species and potential
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();

	const BatterySpeciesVCMap& svcmap = *_svcMapVector[0];
	typename BatterySpeciesVCMap::const_iterator pos = svcmap.find(n);
        if (pos==svcmap.end())
	{   
           throw CException("BatteryModel: Error in Species VC Map");
	}
	const BatterySpeciesVC<T>& svc = *(pos->second);
	const T massDiffusivity = svc["massDiffusivity"]; 

	const TGradArray& potentialGradCell = dynamic_cast<const TGradArray&>(_batteryModelFields.potential_gradient[cells]);
	const TGradArray& speciesGradCell = dynamic_cast<const TGradArray&>(_batteryModelFields.speciesGradient[cells]);
	TArray& heatSource = dynamic_cast<TArray&>(_batteryModelFields.heatSource[cells]);
	for (int c=0; c<cells.getCount(); c++)
	  {
	    const TGradType pGrad =  potentialGradCell[c];
	    const TGradType sGrad =  speciesGradCell[c];
	    T CellSource = pGrad[0]*sGrad[0] + pGrad[1]*sGrad[1];
	    if ((*_meshes[0]).getDimension() == 3)
	      CellSource += pGrad[2]*sGrad[2];
	    //heatSource[c] = CellSource; 
	    heatSource[c] = CellSource*massDiffusivity*96485.0; 
	  }
      } 

    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>(_meshes,_geomFields,
                                            _batteryModelFields.temperature,
                                            _batteryModelFields.thermalConductivity,
                                            _batteryModelFields.temperatureGradient));
    discretizations.push_back(dd);   

    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativeDiscretization<T,T,T>
             (_meshes,_geomFields,
              _batteryModelFields.temperature,
              _batteryModelFields.temperatureN1,
              _batteryModelFields.temperatureN2,
              _batteryModelFields.rhoCp,
              _options["timeStep"]));
        
        discretizations.push_back(td);
    }

    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  _batteryModelFields.temperature,
	  _batteryModelFields.heatSource));
    discretizations.push_back(sd);
 
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
    
    

    

    /* linearize double shell mesh for thermal special interface*/
    
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
	      BatteryLinearizeThermalInterface<T, T, T> lbv (_geomFields,
							_batteryModelFields.temperature,
							sFields.concentration,
							_batteryModelFields.potential,
							_options["ButlerVolmerRRConstant"],
							Anode,
							Cathode);

	      lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	      
	    }
	  else
	    {
	  
	      LinearizeInterfaceJump<T, T, T> lsm (T(1.0),
						   T(0.0),
						   _batteryModelFields.temperature);

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
	    const int nFaces = faces.getCount();
            const BatteryThermalBC<T>& bc = *_tbcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.temperature,
                                  _batteryModelFields.heatFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedTemperature")
            {
            	FloatValEvaluator<T> bT(bc.getVal("specifiedTemperature"), faces);
                for(int f=0; f<nFaces; f++)
		{
                    gbc.applyDirichletBC(f, bT[f]);
		}
            }
            else if (bc.bcType == "SpecifiedHeatFlux")
            {
                const T specifiedFlux(bc["specifiedHeatFlux"]);
                gbc.applyNeumannBC(specifiedFlux);
            }
	    else if (bc.bcType == "Symmetry")
            {
                T zeroFlux(NumTypeTraits<T>::getZero());
                gbc.applyNeumannBC(zeroFlux);
            }
	    else
                throw CException(bc.bcType + " not implemented for Thermal in BatteryModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
	  
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

	    
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.temperature,
                                  _batteryModelFields.heatFlux, 
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }


void linearizePC_Thermal(LinearSystem& ls)
  {

    // only lithium (species 0) for now
    const BatterySpeciesVCMap& svcmap = *_svcMapVector[0];
    const BatterySpeciesBCMap& sbcmap = *_sbcMapVector[0];

    // get combined gradients of potential and species
    GradientModel<VectorT3> pstGradientModel(_meshes,_batteryModelFields.potentialSpeciesTemp,
					  _batteryModelFields.potentialSpeciesTempGradient,_geomFields);
    pstGradientModel.compute();

    //populate lnSpeciesConc fields for inclusion of ln term in potential equation
    bool foundElectrolyte = false;

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	if (n==_options["BatteryElectrolyteMeshID"])
	  {
	    foundElectrolyte = true;
	    const Mesh& mesh = *_meshes[n];
	    const StorageSite& cells = mesh.getCells();
	    const VectorT3Array& speciesPotentialTemp = dynamic_cast<const VectorT3Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    TArray& lnSpeciesConc = dynamic_cast<TArray&>(_batteryModelFields.lnLithiumConcentration[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		T CellConc = (speciesPotentialTemp[c])[1];
		if (CellConc <= 0.0)
		  {
		    cout << "Error: Cell Concentration <= 0   MeshID: " << n << " CellNum: " << c << endl;
		    CellConc = 0.01;
		  }
		lnSpeciesConc[c] = log(CellConc); 
	      }
	  }
      }

    if ((!(foundElectrolyte))&&(_options.ButlerVolmer))
      cout << "Warning: Electrolyte Mesh ID not set." << endl;

    //compute gradient for ln term discretization
    GradientModel<T> lnSpeciesGradientModel(_meshes,_batteryModelFields.lnLithiumConcentration,
					    _batteryModelFields.lnLithiumConcentrationGradient,_geomFields);
    lnSpeciesGradientModel.compute();

    // fill source field with dot product of gradients from species and potential
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();

	typename BatterySpeciesVCMap::const_iterator pos = svcmap.find(n);
	if (pos==svcmap.end())
	  {   
	    throw CException("BatteryModel: Error in Species VC Map");
	  }
	const BatterySpeciesVC<T>& svc = *(pos->second);
	const T massDiffusivity = svc["massDiffusivity"]; 

	const VectorT3GradArray& pstGradCell = dynamic_cast<const VectorT3GradArray&>(_batteryModelFields.potentialSpeciesTempGradient[cells]);
	TArray& heatSource = dynamic_cast<TArray&>(_batteryModelFields.heatSource[cells]);
	for (int c=0; c<cells.getCount(); c++)
	  {
	    const VectorT3Grad combinedCellGradient = pstGradCell[c];
	    TGradType pGrad(NumTypeTraits<TGradType>::getZero());
	    TGradType sGrad(NumTypeTraits<TGradType>::getZero());
	    pGrad[0] = combinedCellGradient[0][0];
	    pGrad[1] = combinedCellGradient[1][0];
	    sGrad[0] = combinedCellGradient[0][1];
	    sGrad[1] = combinedCellGradient[1][1];
	    T CellSource = pGrad[0]*sGrad[0] + pGrad[1]*sGrad[1];
	    if ((*_meshes[0]).getDimension() == 3)
	      {
		pGrad[2] = combinedCellGradient[2][0];
		sGrad[2] = combinedCellGradient[2][1];
		CellSource += pGrad[2]*sGrad[2];
	      }
	    heatSource[c] = CellSource*massDiffusivity*96485.0; 
	  }
      }	   
      
  
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new BatteryPCDiffusionDiscretization<VectorT3,DiagTensorT3,DiagTensorT3>
	 (_realMeshes,_geomFields,
	  _batteryModelFields.potentialSpeciesTemp,
	  _batteryModelFields.potentialSpeciesTempDiffusivity,
	  _batteryModelFields.potentialSpeciesTempGradient));
    discretizations.push_back(dd);

    //discretize ln term (only affects electrolyte mesh of potential model)
    shared_ptr<Discretization>
      bedd(new BatteryPCBinaryElectrolyteDiscretization<VectorT3,DiagTensorT3,DiagTensorT3>(_realMeshes,_geomFields,
							     _batteryModelFields.potentialSpeciesTemp,
							     _batteryModelFields.potentialSpeciesTempDiffusivity,
							     _batteryModelFields.lnLithiumConcentration,
							     _batteryModelFields.lnLithiumConcentrationGradient,
							     _options.thermalModelPC,
							     _options["BatteryElectrolyteMeshID"]));
    discretizations.push_back(bedd); 
    
    if (_options.transient)
      {
        shared_ptr<Discretization>
          td(new BatteryPCTimeDerivativeDiscretization<VectorT3,DiagTensorT3,DiagTensorT3>
             (_realMeshes,_geomFields,
              _batteryModelFields.potentialSpeciesTemp,
              _batteryModelFields.potentialSpeciesTempN1,
              _batteryModelFields.potentialSpeciesTempN2,
	      _batteryModelFields.rhoCp,
	      _options.thermalModelPC,
              _options["timeStep"]));
        
        discretizations.push_back(td);
      }
    
    if (1)
      {
        shared_ptr<Discretization>
          sd(new BatteryPCHeatSourceDiscretization<VectorT3>
             (_realMeshes,_geomFields,
              _batteryModelFields.potentialSpeciesTemp,
	      _batteryModelFields.heatSource));
        
        discretizations.push_back(sd);
      }
    
    Linearizer linearizer;

    linearizer.linearize(discretizations,_realMeshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

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
		
		BatteryPCLinearizeInterface_BV<VectorT3, SquareTensorT3, SquareTensorT3, DiagTensorT3> lbv (_geomFields,
							_options["ButlerVolmerRRConstant"],
							_options["interfaceSpeciesUnderRelax"],		     	      
							Anode,
							Cathode,
							_options.thermalModelPC,			       
							_batteryModelFields.potentialSpeciesTemp);

							lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );

	      }
	    else
	      {
		LinearizeInterfaceJump<VectorT3, SquareTensorT3, SquareTensorT3> lsm (T(1.0),
										      NumTypeTraits<VectorT3>::getZero(),
						     _batteryModelFields.potentialSpeciesTemp);

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

	    const BatteryThermalBC<T>& tbc = *_tbcMap[fg.id];

            BatteryPC_BCS<VectorT3,DiagTensorT3,DiagTensorT3> bbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.potentialSpeciesTemp,
                                  _batteryModelFields.potentialSpeciesTempFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (sbc.bcType == "SpecifiedConcentration")
            {
                FloatValEvaluator<T>
                  sbT(sbc.getVal("specifiedConcentration"),faces);
    
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
	    
	    if (tbc.bcType == "SpecifiedTemperature")
	      {
		FloatValEvaluator<T> tbT(tbc.getVal("specifiedTemperature"),faces);
    
		bbc.applySingleEquationDirichletBC(tbT,2);
	      }
	    else if (tbc.bcType == "SpecifiedHeatFlux")
	      {
		const T specifiedFlux(tbc["specifiedHeatFlux"]);
		bbc.applySingleEquationNeumannBC(specifiedFlux,2);
	      }
	    else if ((tbc.bcType == "Symmetry"))
	      {
		T zeroFlux(NumTypeTraits<T>::getZero());
		bbc.applySingleEquationNeumannBC(zeroFlux,2);
	      }
	    else
	      throw CException(tbc.bcType + " not implemented for temperature in BatteryModel");
	      
        }
	
	if (mesh.isDoubleShell())
	  {
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
		GenericBCS<VectorT3,SquareTensorT3,SquareTensorT3> gbc(faces,mesh,
								       _geomFields,
								       _batteryModelFields.potentialSpeciesTemp,
								       _batteryModelFields.potentialSpeciesTempFlux,
								       ls.getMatrix(), ls.getX(), ls.getB());

		gbc.applyInterfaceBC();
	      }
	  }
	else
	  {
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
		GenericBCS<VectorT3,DiagTensorT3,DiagTensorT3> gbc(faces,mesh,
								       _geomFields,
								       _batteryModelFields.potentialSpeciesTemp,
								       _batteryModelFields.potentialSpeciesTempFlux,
								       ls.getMatrix(), ls.getX(), ls.getB());

		gbc.applyInterfaceBC();
	      }
	  }
	//set all temperature residuals to zero if no thermal model PC (comes in to play if trying to still do advanceThermal separatly)
/*
	if (!_options.thermalModelPC)
	  {
	    for (int n=0; n<numMeshes; n++)
	      {
		const Mesh& mesh = *_meshes[n];
		const StorageSite& cells = mesh.getCells();
		const MultiField::ArrayIndex cVarIndex(&_batteryModelFields.potentialSpeciesTemp,&cells);
		VectorT3Array& rCell = dynamic_cast<VectorT3Array&>((ls.getB())[cVarIndex]);
		for (int c=0; c<cells.getCount(); c++)
		  {
		    (rCell[c])[2] = 0.0;
		  }
	      }
	  }*/
    }
  }

void linearizePC(LinearSystem& ls)
  {
    // only lithium (species 0) for now
    const BatterySpeciesBCMap& sbcmap = *_sbcMapVector[0];

    // get combined gradients of potential and species
    GradientModel<VectorT2> pstGradientModel(_meshes,_batteryModelFields.potentialSpeciesTemp,
					  _batteryModelFields.potentialSpeciesTempGradient,_geomFields);
    pstGradientModel.compute();

    //populate lnSpeciesConc fields for inclusion of ln term in potential equation
    bool foundElectrolyte = false;

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	if (n==_options["BatteryElectrolyteMeshID"])
	  {
	    foundElectrolyte = true;
	    const Mesh& mesh = *_meshes[n];
	    const StorageSite& cells = mesh.getCells();
	    const VectorT2Array& speciesPotentialTemp = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    TArray& lnSpeciesConc = dynamic_cast<TArray&>(_batteryModelFields.lnLithiumConcentration[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		T CellConc = (speciesPotentialTemp[c])[1];
		if (CellConc <= 0.0)
		  {
		    cout << "Error: Cell Concentration <= 0   MeshID: " << n << " CellNum: " << c << endl;
		    CellConc = 0.01;
		  }
		lnSpeciesConc[c] = log(CellConc); 
	      }
	  }
      }

    if ((!(foundElectrolyte))&&(_options.ButlerVolmer))
      cout << "Warning: Electrolyte Mesh ID not set." << endl;

    //compute gradient for ln term discretization
    GradientModel<T> lnSpeciesGradientModel(_meshes,_batteryModelFields.lnLithiumConcentration,
					    _batteryModelFields.lnLithiumConcentrationGradient,_geomFields);
    lnSpeciesGradientModel.compute();

    
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new BatteryPCDiffusionDiscretization<VectorT2,DiagTensorT2,DiagTensorT2>
	 (_realMeshes,_geomFields,
	  _batteryModelFields.potentialSpeciesTemp,
	  _batteryModelFields.potentialSpeciesTempDiffusivity,
	  _batteryModelFields.potentialSpeciesTempGradient));
    discretizations.push_back(dd);

    //discretize ln term (only affects electrolyte mesh of potential model)
    shared_ptr<Discretization>
      bedd(new BatteryPCBinaryElectrolyteDiscretization<VectorT2,DiagTensorT2,DiagTensorT2>(_realMeshes,_geomFields,
							     _batteryModelFields.potentialSpeciesTemp,
							     _batteryModelFields.potentialSpeciesTempDiffusivity,
							     _batteryModelFields.lnLithiumConcentration,
							     _batteryModelFields.lnLithiumConcentrationGradient,
							     _options.thermalModelPC,
							     _options["BatteryElectrolyteMeshID"]));
    discretizations.push_back(bedd); 
    
    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new BatteryPCTimeDerivativeDiscretization<VectorT2,DiagTensorT2,DiagTensorT2>
             (_realMeshes,_geomFields,
              _batteryModelFields.potentialSpeciesTemp,
              _batteryModelFields.potentialSpeciesTempN1,
              _batteryModelFields.potentialSpeciesTempN2,
	      _batteryModelFields.rhoCp,
	      _options.thermalModelPC,
              _options["timeStep"]));
        
        discretizations.push_back(td);
	}
    
    Linearizer linearizer;

    linearizer.linearize(discretizations,_realMeshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

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
		
		BatteryPCLinearizeInterface_BV<VectorT2, SquareTensorT2, SquareTensorT2, DiagTensorT2> lbv (_geomFields,
							_options["ButlerVolmerRRConstant"],
							_options["interfaceSpeciesUnderRelax"],		     	      
							Anode,
							Cathode,
							_options.thermalModelPC,			       
							_batteryModelFields.potentialSpeciesTemp);

							lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );

	      }
	    else
	      {
		LinearizeInterfaceJump<VectorT2, SquareTensorT2, SquareTensorT2> lsm (T(1.0),
										      NumTypeTraits<VectorT2>::getZero(),
						     _batteryModelFields.potentialSpeciesTemp);

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

            BatteryPC_BCS<VectorT2,DiagTensorT2,DiagTensorT2> bbc(faces,mesh,
                                  _geomFields,
                                  _batteryModelFields.potentialSpeciesTemp,
                                  _batteryModelFields.potentialSpeciesTempFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (sbc.bcType == "SpecifiedConcentration")
            {
                FloatValEvaluator<T>
                  sbT(sbc.getVal("specifiedConcentration"),faces);
    
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
	
	if (mesh.isDoubleShell())
	  {
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
		GenericBCS<VectorT2,SquareTensorT2,SquareTensorT2> gbc(faces,mesh,
								       _geomFields,
								       _batteryModelFields.potentialSpeciesTemp,
								       _batteryModelFields.potentialSpeciesTempFlux,
								       ls.getMatrix(), ls.getX(), ls.getB());

		gbc.applyInterfaceBC();
	      }
	  }
	else
	  {
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	      {
		const FaceGroup& fg = *fgPtr;
		const StorageSite& faces = fg.site;
		GenericBCS<VectorT2,DiagTensorT2,DiagTensorT2> gbc(faces,mesh,
								       _geomFields,
								       _batteryModelFields.potentialSpeciesTemp,
								       _batteryModelFields.potentialSpeciesTempFlux,
								       ls.getMatrix(), ls.getX(), ls.getB());

		gbc.applyInterfaceBC();
	      }
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
              dynamic_cast<const TArray&>(_batteryModelFields.heatFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += heatFlux[f];
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
	  const TArray& heatFlux =
	    dynamic_cast<const TArray&>(_batteryModelFields.heatFlux[faces]);
	  for(int f=0; f<nFaces; f++)
	    r += heatFlux[f];
	  found=true;
        }
    }


    if (!found)
      throw CException("getPotentialFluxIntegral: invalid faceGroupID");
    return r;
  }

T getAverageConcentration(const Mesh& mesh, const int m)
  {
    BatterySpeciesFields& sFields = *_speciesFieldsVector[m];
    const StorageSite& cells = mesh.getCells();
    const TArray& concCell = dynamic_cast<const TArray&>(sFields.concentration[cells]);
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    T totalVolume = 0.0;
    T weightedConcentration = 0.0;
    for(int c=0; c<cells.getSelfCount(); c++)
      {
	totalVolume += cellVolume[c];
	weightedConcentration += concCell[c]*cellVolume[c];
      }
    return weightedConcentration/totalVolume;
  }


 void printMatrixElementsOnFace(const Mesh& mesh, const int fgId, LinearSystem& ls, Field& varField)
 {
   
   const FaceGroup& fg = mesh.getFaceGroup(fgId);
   const StorageSite& faces = fg.site;
   const StorageSite& cells = mesh.getCells();
   const int nFaces = faces.getCount();
   const CRConnectivity& faceCells = mesh.getFaceCells(faces);

   MultiFieldMatrix& mfmatrix = ls.getMatrix();
   MultiField& rField = ls.getB();
   MultiField& xField = ls.getX();

   const MultiField::ArrayIndex cVarIndex(&varField,&cells);
   CRMatrix<T,T,T>& matrix = dynamic_cast<CRMatrix<T,T,T>&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
   TArray& rCell = dynamic_cast<TArray&>(rField[cVarIndex]);
   TArray& xCell = dynamic_cast<TArray&>(xField[cVarIndex]);
   typename CRMatrix<T,T,T>::DiagArray& diag = matrix.getDiag();

   for(int f=0; f<nFaces; f++)
     {
       const int c0 = faceCells(f,0);
       const int c1 = faceCells(f,1);
       T& offdiagC0_C1 = matrix.getCoeff(c0,  c1);
       T& offdiagC1_C0 = matrix.getCoeff(c1,  c0);
       const T rC0= rCell[c0];
       const T rC1= rCell[c1];
       const T diagC0= diag[c0];
       const T diagC1= diag[c1];
       const T xC0 = xCell[c0];
       const T xC1 = xCell[c1];

       if (f==0)
	 {
	   cout << "--------------" << endl;
	   cout << mesh.getID() << ":" << f << endl;
	   cout << offdiagC0_C1 << endl;
	   cout << offdiagC1_C0 << endl;
	   cout << rC0 << endl;
	   cout << rC1 << endl;
	   cout << diagC0 << endl;
	   cout << diagC1 << endl;
	   cout << xC0 << endl;
	   cout << xC1 << endl;
	   
	   cout << "--------------" << endl;

	 }
       
     }
   return;
 }

 LinearSystem& matrixMerge(LinearSystem& ls, const int m)
 {
   /*
   typedef CRMatrix<double,double,double>  CombinedMatrix;
   const MultiFieldMatrix& mfMatrix = ls.getMatrix();
   const MultiField& bField = ls.getB();
   const MultiField& deltaField = ls.getDelta();
   if (bField.getLength() == 1)
      return ls;

    const int numMeshes = _meshes.size();
    int totalBoundaryCells = 0;
    int totalInteriorCells = 0;
    for (int n1=0; n1<numMeshes; n1++)
    {
        Mesh& mesh = *_meshes[n1];
	StorageSite& cells = mesh.getCells();
        //const CRConnectivity& cellCells = mesh.getCellCells();
	
	
	int n1InterfaceCells = 0;
	for (int n2=0; n2<numMeshes; n2++)
	  {
	    StorageSite& otherCells = (*_meshes[n2]).getCells();
	    StorageSite::ScatterMap& n1ScatterMap = cells.getScatterMap();

	    typename StorageSite::ScatterMap::const_iterator pos = n1ScatterMap.find(&otherCells);
            if (pos==n1ScatterMap.end())
	      {   
		//no interfaces with mesh n2
	      }
	    else
	      {
		shared_ptr< Array<int> > n1n2ScatterOrigPtr = n1ScatterMap[&otherCells];
		n1InterfaceCells += n1n2ScatterOrigPtr->getLength();
	      }

	  //shared_ptr<CRConnectivity> flatConnPtr = origConn.getMultiTranspose(blockSize);

	    //CombinedMatrix combinedMatrix(*combinedConnPtr);

	    // origMatrix.setFlatMatrix(flatMatrix);
	    //ls.getMatrix().addMatrix(tIndex,tIndex,m);
	  }
	const int n1InteriorCells = cells.getSelfCount();
	const int n1BoundaryCells = cells.getCount() - n1InterfaceCells - n1InteriorCells;
	totalBoundaryCells += n1BoundaryCells;
	totalInteriorCells += n1InteriorCells;
	}

    cout << totalBoundaryCells << endl;
    cout << totalInteriorCells << endl;

    const MultiField::ArrayIndex rowIndex = bField.getArrayIndex(0);
    const Matrix& origMatrix = mfMatrix.getMatrix(rowIndex,rowIndex);
    const CRConnectivity& origConn = origMatrix.getConnectivity();
    const ArrayBase& origB = bField[rowIndex];

    LinearSystem lsCombined;
   
    shared_ptr<CRConnectivity> trPtr(new CRConnectivity(origConn.getColSite(),origConn.getRowSite()));
    CRConnectivity& tr = *trPtr;

    tr.getrowDim *= varSize;
    tr._colDim *= varSize;
  
    tr.initCount();

    const Array<int>& myRow = *_row;
    const Array<int>& myCol = *_col;*/
   

   return ls; 
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

        LinearSystem ls,lsShell;
        //LinearSystem ls;

	initSpeciesLinearization(ls, lsShell, m);
        //initSpeciesLinearization(ls, ls, m);
        
        ls.initAssembly();
        lsShell.initAssembly();

	linearizeSpecies(ls, lsShell, m);
        //linearizeSpecies(ls, ls, m);

	ls.initSolve();

	//cout << "BEFORE LS" << endl;
	///BatterySpeciesFields& sFields = *_speciesFieldsVector[m];
	//printMatrixElementsOnFace(*_meshes[0], 8, ls, sFields.concentration);
	//printMatrixElementsOnFace(*_meshes[1], 11, ls, sFields.concentration);


        MFRPtr rNorm(_options.getLinearSolverSpecies().solve(ls));

        if (!iNorm) iNorm = rNorm;        
        MFRPtr normRatio((*rNorm)/(*iNorm));

        //cout << "Species Number: " << m << endl;
	if (((_niters+1) % _options.advanceVerbosity)==0)
	  cout << _niters << ": " << *rNorm << endl;
        
        _options.getLinearSolverSpecies().cleanup();

        ls.postSolve();
        ls.updateSolution();
	
	//cout << "AFTER LS" << endl;
	//printMatrixElementsOnFace(*_meshes[0], 8, ls, sFields.concentration);
	//printMatrixElementsOnFace(*_meshes[1], 11, ls, sFields.concentration);
	

	//updateShellGhosts();  not useful right now
	lsShell.initSolve();
	MFRPtr rNorm2(_options.getLinearSolverSpecies().solve(lsShell));
        _options.getLinearSolverSpecies().cleanup();
	lsShell.postSolve();
	lsShell.updateSolution();
	cout << _niters << ": " << *rNorm2 << endl;
	


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

	if ((_niters % _options.advanceVerbosity)==0)
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

void advanceThermal(const int niter)
  {
    for(int n=0; n<niter; n++)
    { 
        LinearSystem ls;

        initThermalLinearization(ls);
        
        ls.initAssembly();

        linearizeThermal(ls);

        ls.initSolve();

        MFRPtr rNorm(_options.getLinearSolverThermal().solve(ls));

        if (!_initialThermalNorm) _initialThermalNorm = rNorm;
        
	_currentThermalResidual = rNorm;

        MFRPtr normRatio((*rNorm)/(*_initialThermalNorm));

	if ((_niters % _options.advanceVerbosity)==0)
	  cout << "Thermal: " << _niters << ": " << *rNorm << endl;
        
        _options.getLinearSolverThermal().cleanup();

        ls.postSolve();
        ls.updateSolution();

        _niters++;

        if (*rNorm < _options.absoluteThermalTolerance ||
            *normRatio < _options.relativeThermalTolerance)
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

       if (_options.thermalModelPC)
	 linearizePC_Thermal(ls);
       else
	 linearizePC(ls);

       ls.initSolve();

       MFRPtr rNorm = _options.getLinearSolverPC().solve(ls);

       if (!_initialPCNorm) _initialPCNorm = rNorm;

       _currentPCResidual = rNorm;

       MFRPtr normRatio((*rNorm)/(*_initialPCNorm));

       _options.getLinearSolverPC().cleanup();

       ls.postSolve();
    
       ls.updateSolution();
       
       if ((_niters % _options.advanceVerbosity)==0)
	 cout << "Point-Coupled: " << _niters << ": " << *rNorm << endl;

       _niters++;

       if (*rNorm < _options.absolutePCTolerance ||
	   *normRatio < _options.relativePCTolerance)
	 break;
     }
 }

  T getSpeciesResidual(const int speciesId)
  {
    MFRPtr& rCurrent = *_currentSpeciesResidual[speciesId];
    BatterySpeciesFields& sFields = *_speciesFieldsVector[speciesId];
    const Field& concentrationField = sFields.concentration;
    ArrayBase& residualArray = (*rCurrent)[concentrationField];
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

  T getThermalResidual()
  {
    const Field& thermalField = _batteryModelFields.temperature;
    ArrayBase& residualArray = (*_currentThermalResidual)[thermalField];
    T *arrayData = (T*)(residualArray.getData());
    const T residual = arrayData[0];
    return residual;
  }

  T getPCResidual(const int v)
  {
    const Field& potentialSpeciesTempField = _batteryModelFields.potentialSpeciesTemp;
    ArrayBase& residualArray = (*_currentPCResidual)[potentialSpeciesTempField];
    T *arrayData = (T*)(residualArray.getData());
    if (v<2)
      {
	const T residual = arrayData[v];
	return residual;
      }
    else if ((v==2)&&(_options.thermalModelPC))
      {
	const T residual = arrayData[v];
	return residual;
      }
    else
      {
	const T residual = 0.0;
	return residual;
      }
  }
  
  void updateShellGhosts()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	if ((mesh.isDoubleShell())&&(!(mesh.isConnectedShell())))
	  {
	    BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium

	    const StorageSite& cells = mesh.getCells();
	    const StorageSite& faces = mesh.getParentFaceGroupSite();
	    const CRConnectivity& cellCells = mesh.getCellCells();
	    TArray& speciesCell = dynamic_cast<TArray&>(sFields.concentration[cells]);

	    const int parentMeshID = mesh.getParentMeshID();
	    const int otherMeshID = mesh.getOtherMeshID();
	    const Mesh& parentMesh = *_meshes[parentMeshID];
	    const Mesh& otherMesh = *_meshes[otherMeshID];

	    const CRConnectivity& parentFaceCells = parentMesh.getFaceCells(faces);
	    const StorageSite& parentCells = parentMesh.getCells();
	    const StorageSite& otherFaces = mesh.getOtherFaceGroupSite();
	    const CRConnectivity& otherFaceCells = otherMesh.getFaceCells(otherFaces);
	    const StorageSite& otherCells = otherMesh.getCells();
	    const TArray& speciesCellParent = dynamic_cast< const TArray&>(sFields.concentration[parentCells]);
	    const TArray& speciesCellOther = dynamic_cast< const TArray&>(sFields.concentration[otherCells]);


	    for (int f=0; f<faces.getCount(); f++)
	      {
		int c0p = parentFaceCells(f,0);
		int c1p = parentFaceCells(f,1);
		if (c1p < parentCells.getSelfCount())
		  { 		   
		    int temp = c0p;
		    c0p = c1p;
		    c1p = temp;
		  }

		int c0o = otherFaceCells(f,0);
		int c1o = otherFaceCells(f,1);
		if (c1o < otherCells.getSelfCount())
		  { 
		    int temp = c0o;
		    c0o = c1o;
		    c1o = temp;
		  }

		const int c0 = f;
		const int c1 = cellCells(f,0);
		const int c2 = cellCells(f,1);
		const int c3 = cellCells(f,2);

		if (c0==0)
		  {
		    cout << speciesCell[c0] << " " << speciesCell[c1] << " " << speciesCell[c2] << " " << speciesCell[c3] << endl;
		  }

		// copy sln varialbe values from interior cells of meshes to ghost cells of shell mesh
		speciesCell[c3] = speciesCellParent[c0p];
		speciesCell[c2] = speciesCellOther[c0o];

		if (c0==0)
		  {
		    cout << speciesCell[c0] << " " << speciesCell[c1] << " " << speciesCell[c2] << " " << speciesCell[c3] << endl;
		  }
	      }
	  }
      }
  }

  void copySeparateToCoupled()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	//copy concentration and potential and thermal
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium
	const TArray& mFSeparate = dynamic_cast<const TArray&>(sFields.concentration[cells]);
	const TArray& pSeparate = dynamic_cast<const TArray&>(_batteryModelFields.potential[cells]);
	const TArray& tSeparate = dynamic_cast<const TArray&>(_batteryModelFields.temperature[cells]);

	if (_options.thermalModelPC)
	  {
	    VectorT3Array& coupledValues = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		T concentration = mFSeparate[c];
		T potential = pSeparate[c];
		T temperature = tSeparate[c];
		VectorT3& coupledCellVector = coupledValues[c];

		//populate coupled field with separate field data
		coupledCellVector[0] = potential;
		coupledCellVector[1] = concentration;
		coupledCellVector[2] = temperature;
	      }
	  }
	else
	  {
	    VectorT2Array& coupledValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		T concentration = mFSeparate[c];
		T potential = pSeparate[c];
		VectorT2& coupledCellVector = coupledValues[c];

		//populate coupled field with separate field data
		coupledCellVector[0] = potential;
		coupledCellVector[1] = concentration;
	      }
	  }

	//copy fluxes
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const TArray& mFFluxSeparate = dynamic_cast<const TArray&>(sFields.massFlux[faces]);
	    const TArray& pFluxSeparate = dynamic_cast<const TArray&>(_batteryModelFields.potential_flux[faces]);  
	    const TArray& tFluxSeparate = dynamic_cast<const TArray&>(_batteryModelFields.heatFlux[faces]);
	    if (_options.thermalModelPC)
	      {
		VectorT3Array& coupledFluxValues = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    T concentrationFlux = mFFluxSeparate[f];
		    T potentialFlux = pFluxSeparate[f];
		    T heatFlux = tFluxSeparate[f];
		    VectorT3& coupledCellFluxVector = coupledFluxValues[f];

		    //populate coupled field with separate field data
		    coupledCellFluxVector[0] = potentialFlux;
		    coupledCellFluxVector[1] = concentrationFlux;
		    coupledCellFluxVector[2] = heatFlux;		      
		  }	
	      }      
	    else
	      {
		VectorT2Array& coupledFluxValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    T concentrationFlux = mFFluxSeparate[f];
		    T potentialFlux = pFluxSeparate[f];
		    VectorT2& coupledCellFluxVector = coupledFluxValues[f];

		    //populate coupled field with separate field data
		    coupledCellFluxVector[0] = potentialFlux;
		    coupledCellFluxVector[1] = concentrationFlux;
		  }	
	      }      
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const TArray& mFFluxSeparate = dynamic_cast<const TArray&>(sFields.massFlux[faces]);
	    const TArray& pFluxSeparate = dynamic_cast<const TArray&>(_batteryModelFields.potential_flux[faces]);  
	    const TArray& tFluxSeparate = dynamic_cast<const TArray&>(_batteryModelFields.heatFlux[faces]);  
	    if (_options.thermalModelPC)
	      {
		VectorT3Array& coupledFluxValues = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    T concentrationFlux = mFFluxSeparate[f];
		    T potentialFlux = pFluxSeparate[f];
		    T heatFlux = tFluxSeparate[f];
		    VectorT3& coupledCellFluxVector = coupledFluxValues[f];

		    //populate coupled field with separate field data
		    coupledCellFluxVector[0] = potentialFlux;
		    coupledCellFluxVector[1] = concentrationFlux;
		    coupledCellFluxVector[2] = heatFlux;	     
		  }	
	      }      
	    else
	      {
		VectorT2Array& coupledFluxValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    T concentrationFlux = mFFluxSeparate[f];
		    T potentialFlux = pFluxSeparate[f];
		    VectorT2& coupledCellFluxVector = coupledFluxValues[f];

		    //populate coupled field with separate field data
		    coupledCellFluxVector[0] = potentialFlux;
		    coupledCellFluxVector[1] = concentrationFlux;
		  }	
	      }    
	  }

	// copy to N1 and N2 if transient
	if (_options.transient)
	  {
	    const TArray& potSeparateN1 = dynamic_cast<const TArray&>(_batteryModelFields.potentialN1[cells]);
	    const TArray& mFSeparateN1 = dynamic_cast<const TArray&>(sFields.concentrationN1[cells]);
	    const TArray& tSeparateN1 = dynamic_cast<const TArray&>(_batteryModelFields.temperatureN1[cells]);
	    if (_options.thermalModelPC)
	      {
		VectorT3Array& coupledValuesN1 = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);

		for (int c=0; c<cells.getCount(); c++)
		  {
		    const T potentialN1 = potSeparateN1[c];
		    const T concentrationN1 = mFSeparateN1[c];
		    const T temperatureN1 = tSeparateN1[c];
		    VectorT3& coupledCellVectorN1 = coupledValuesN1[c];

		    //populate coupled field with separate field data
		    coupledCellVectorN1[1] = concentrationN1;
		    coupledCellVectorN1[0] = potentialN1;
		    coupledCellVectorN1[2] = temperatureN1;
		  }

		if (_options.timeDiscretizationOrder > 1)
		  {
		    VectorT3Array& coupledValuesN2 = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		    const TArray& mFSeparateN2 = dynamic_cast<const TArray&>(sFields.concentrationN2[cells]);
		    const TArray& tSeparateN2 = dynamic_cast<const TArray&>(_batteryModelFields.temperatureN2[cells]);
		    for (int c=0; c<cells.getCount(); c++)
		      {
			const T concentrationN2 = mFSeparateN2[c];
			const T temperatureN2 = tSeparateN2[c];
			VectorT3& coupledCellVectorN2 = coupledValuesN2[c];

			//populate coupled field with separate field data
			coupledCellVectorN2[1] = concentrationN2;
			coupledCellVectorN2[2] = temperatureN2;	 
		      }
		  }
	      }
	    else
	      {
		VectorT2Array& coupledValuesN1 = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);

		for (int c=0; c<cells.getCount(); c++)
		  {
		    const T potentialN1 = potSeparateN1[c];
		    const T concentrationN1 = mFSeparateN1[c];
		    VectorT2& coupledCellVectorN1 = coupledValuesN1[c];

		    //populate coupled field with separate field data
		    coupledCellVectorN1[1] = concentrationN1;
		    coupledCellVectorN1[0] = potentialN1;
		  }

		if (_options.timeDiscretizationOrder > 1)
		  {
		    VectorT2Array& coupledValuesN2 = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		    const TArray& mFSeparateN2 = dynamic_cast<const TArray&>(sFields.concentrationN2[cells]);
		    for (int c=0; c<cells.getCount(); c++)
		      {
			const T concentrationN2 = mFSeparateN2[c];
			VectorT2& coupledCellVectorN2 = coupledValuesN2[c];

			//populate coupled field with separate field data
			coupledCellVectorN2[1] = concentrationN2;
		      }
		  }
	      }
	  }
      }
  }

  void copyCoupledToSeparate()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium

	// copy conc and potential
	TArray& mFSeparate = dynamic_cast<TArray&>(sFields.concentration[cells]);
	TArray& pSeparate = dynamic_cast<TArray&>(_batteryModelFields.potential[cells]);
	TArray& tSeparate = dynamic_cast<TArray&>(_batteryModelFields.temperature[cells]);
	if (_options.thermalModelPC)
	  {
	    const VectorT3Array& coupledValues = dynamic_cast<const VectorT3Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		VectorT3 coupledCellVector = coupledValues[c];
		T potential = coupledCellVector[0];
		T concentration = coupledCellVector[1];
		T temperature = coupledCellVector[2];

		//populate separate fields with coupled field data
		pSeparate[c] = potential;
		mFSeparate[c] = concentration;
		tSeparate[c] = temperature;
	      }
	  }
	else
	  {
	    const VectorT2Array& coupledValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialSpeciesTemp[cells]);
	    for (int c=0; c<cells.getCount(); c++)
	      {
		VectorT2 coupledCellVector = coupledValues[c];
		T potential = coupledCellVector[0];
		T concentration = coupledCellVector[1];
		
		//populate separate fields with coupled field data
		pSeparate[c] = potential;
		mFSeparate[c] = concentration;
	      }
	  }
	// copy fluxes
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    TArray& mFFluxSeparate = dynamic_cast<TArray&>(sFields.massFlux[faces]);
	    TArray& pFluxSeparate = dynamic_cast<TArray&>(_batteryModelFields.potential_flux[faces]); 
	    TArray& tFluxSeparate = dynamic_cast<TArray&>(_batteryModelFields.heatFlux[faces]);
	    if (_options.thermalModelPC)
	      {
		const VectorT3Array& coupledFluxValues = dynamic_cast<const VectorT3Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    VectorT3 coupledCellFluxVector = coupledFluxValues[f];
		    T potentialFlux = coupledCellFluxVector[0];
		    T concentrationFlux = coupledCellFluxVector[1];
		    T heatFlux = coupledCellFluxVector[2];	     

		    //populate separate fields with coupled field data
		    pFluxSeparate[f] = potentialFlux;
		    mFFluxSeparate[f] = concentrationFlux;
		    tFluxSeparate[f] = heatFlux;
		  }	 
	      }     
	    else
	      {
		const VectorT2Array& coupledFluxValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    VectorT2 coupledCellFluxVector = coupledFluxValues[f];
		    T potentialFlux = coupledCellFluxVector[0];
		    T concentrationFlux = coupledCellFluxVector[1];

		    //populate separate fields with coupled field data
		    pFluxSeparate[f] = potentialFlux;
		    mFFluxSeparate[f] = concentrationFlux;
		  }	 
	      } 
	  }
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    TArray& mFFluxSeparate = dynamic_cast<TArray&>(sFields.massFlux[faces]);
	    TArray& pFluxSeparate = dynamic_cast<TArray&>(_batteryModelFields.potential_flux[faces]); 
	    TArray& tFluxSeparate = dynamic_cast<TArray&>(_batteryModelFields.heatFlux[faces]);
	    if (_options.thermalModelPC)
	      {
		const VectorT3Array& coupledFluxValues = dynamic_cast<const VectorT3Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    VectorT3 coupledCellFluxVector = coupledFluxValues[f];
		    T potentialFlux = coupledCellFluxVector[0];
		    T concentrationFlux = coupledCellFluxVector[1];
		    T heatFlux = coupledCellFluxVector[2];
		     
		    //populate separate fields with coupled field data
		    pFluxSeparate[f] = potentialFlux;
		    mFFluxSeparate[f] = concentrationFlux;
		    tFluxSeparate[f] = heatFlux;
		      
		  }	
	      }    
	    else
	      {
		const VectorT2Array& coupledFluxValues = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialSpeciesTempFlux[faces]);
		for (int f=0; f<faces.getCount(); f++)
		  {
		    VectorT2 coupledCellFluxVector = coupledFluxValues[f];
		    T potentialFlux = coupledCellFluxVector[0];
		    T concentrationFlux = coupledCellFluxVector[1];
		    
		    //populate separate fields with coupled field data
		    pFluxSeparate[f] = potentialFlux;
		    mFFluxSeparate[f] = concentrationFlux;
		  }	
	      }    
	  }

	// copy to N1 and N2 if transient
	if (_options.transient)
	  {
	    TArray& potSeparateN1 = dynamic_cast<TArray&>(_batteryModelFields.potentialN1[cells]);
	    TArray& mFSeparateN1 = dynamic_cast<TArray&>(sFields.concentrationN1[cells]);
	    TArray& tSeparateN1 = dynamic_cast<TArray&>(_batteryModelFields.temperatureN1[cells]);

	    if (_options.thermalModelPC)
	      {
		const VectorT3Array& coupledValuesN1 = dynamic_cast< const VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);

		for (int c=0; c<cells.getCount(); c++)
		  {
		    const VectorT3& coupledCellVectorN1 = coupledValuesN1[c];

		    //populate separate field with coupled field data
		    mFSeparateN1[c] = coupledCellVectorN1[1];
		    potSeparateN1[c] = coupledCellVectorN1[0];
		    tSeparateN1[c] = coupledCellVectorN1[2]; 
		  }

		if (_options.timeDiscretizationOrder > 1)
		  {
		    const VectorT3Array& coupledValuesN2 = dynamic_cast<const VectorT3Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		    TArray& mFSeparateN2 = dynamic_cast<TArray&>(sFields.concentrationN2[cells]);
		    TArray& tSeparateN2 = dynamic_cast<TArray&>(_batteryModelFields.temperatureN2[cells]);
		    for (int c=0; c<cells.getCount(); c++)
		      {
			const VectorT3& coupledCellVectorN2 = coupledValuesN2[c];

			//populate separate field with coupled field data
			mFSeparateN2[c] = coupledCellVectorN2[1];
			tSeparateN2[c] = coupledCellVectorN2[2];	  
		      }
		  }
	      }
	    else
	      {
		const VectorT2Array& coupledValuesN1 = dynamic_cast< const VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN1[cells]);

		for (int c=0; c<cells.getCount(); c++)
		  {
		    const VectorT2& coupledCellVectorN1 = coupledValuesN1[c];

		    //populate separate field with coupled field data
		    mFSeparateN1[c] = coupledCellVectorN1[1];
		    potSeparateN1[c] = coupledCellVectorN1[0];
		  }

		if (_options.timeDiscretizationOrder > 1)
		  {
		    const VectorT2Array& coupledValuesN2 = dynamic_cast<const VectorT2Array&>(_batteryModelFields.potentialSpeciesTempN2[cells]);
		    TArray& mFSeparateN2 = dynamic_cast<TArray&>(sFields.concentrationN2[cells]);
		    for (int c=0; c<cells.getCount(); c++)
		      {
			const VectorT2& coupledCellVectorN2 = coupledValuesN2[c];

			//populate separate field with coupled field data
			mFSeparateN2[c] = coupledCellVectorN2[1];
		      }
		  }
	      }
	  }
      }
  }

void copyPCDiffusivity()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	//copy concentration diff and potential cond and heat thermal cond
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	BatterySpeciesFields& sFields = *_speciesFieldsVector[0]; //first species is lithium
	const TArray& mFDiffSeparate = dynamic_cast<const TArray&>(sFields.diffusivity[cells]);
	const TArray& pDiffSeparate = dynamic_cast<const TArray&>(_batteryModelFields.conductivity[cells]);
	const TArray& tDiffSeparate = dynamic_cast<const TArray&>(_batteryModelFields.thermalConductivity[cells]);

	if (_options.thermalModelPC)
	  {
	    VectorT3Array& coupledDiffValues = dynamic_cast<VectorT3Array&>(_batteryModelFields.potentialSpeciesTempDiffusivity[cells]);

	    for (int c=0; c<cells.getCount(); c++)
	      {
		T massDiff = mFDiffSeparate[c];
		T potentialDiff = pDiffSeparate[c];
		T thermalDiff = tDiffSeparate[c];
		VectorT3& coupledCellDiffVector = coupledDiffValues[c];

		//populate coupled field with separate field data
		coupledCellDiffVector[0] = potentialDiff;
		coupledCellDiffVector[1] = massDiff;
		coupledCellDiffVector[2] = thermalDiff;
	      }
	  }
	else
	  {
	    VectorT2Array& coupledDiffValues = dynamic_cast<VectorT2Array&>(_batteryModelFields.potentialSpeciesTempDiffusivity[cells]);

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
  }

T getFaceGroupArea(const Mesh& mesh, const int fgID)
 {
   foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
     {
       const FaceGroup& fg = *fgPtr;
       if (fg.id == fgID)
	 {
	   const StorageSite& faces = fg.site;
	   const int nFaces = faces.getCount();
	   const TArray& faceAreaMag = dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	   T totalArea = 0.0;
	   for(int f=0; f<nFaces; f++)
	     {
	       totalArea += faceAreaMag[f];
	       if (f%1000==0)
		 cout << "Percent Done: " << (f*100.0/nFaces) << "%" << endl;
	     }
	   return totalArea;
	 }
     }
   throw CException("getFaceGroupArea: No face group with given id");
 }

T getFaceGroupVoltage(const Mesh& mesh, const int fgID)
 {
   foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
     {
       const FaceGroup& fg = *fgPtr;
       if (fg.id == fgID)
	 {
	   const StorageSite& faces = fg.site;
	   const int nFaces = faces.getCount();
	   const TArray& faceAreaMag = dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	   const StorageSite& cells = mesh.getCells();
	   const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	   const TArray& cellPotential = dynamic_cast<const TArray&>(_batteryModelFields.potential[cells]);
	   T totalArea = 0.0;
	   T voltageSum = 0.0;
	   for(int f=0; f<nFaces; f++)
	     {
	       //const T faceVoltage = harmonicAverage(cellPotential[faceCells(f,0)],cellPotential[faceCells(f,1)]);
	       const T faceVoltage = cellPotential[faceCells(f,1)];
	       totalArea += faceAreaMag[f];
	       voltageSum += faceVoltage*faceAreaMag[f];
	     }
	   return (voltageSum/totalArea);
	 }
     }
   throw CException("getFaceGroupVoltage: No face group with given id");
 }

T getMeshVolume(const Mesh& mesh)
  {
    const StorageSite& cells = mesh.getCells();  
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    T totalVolume = 0.0;
    for(int c=0; c<cells.getSelfCount(); c++)
      {
	totalVolume += cellVolume[c];
      }
    return totalVolume;
  }

private:
  const MeshList _realMeshes;
  const MeshList _meshes;
  GeomFields& _geomFields;
  
  vector<BatterySpeciesFields*> _speciesFieldsVector;

  vector<BatterySpeciesBCMap*> _sbcMapVector;
  vector<BatterySpeciesVCMap*> _svcMapVector;
  BatteryPotentialBCMap _pbcMap;
  BatteryPotentialVCMap _pvcMap;
  BatteryThermalBCMap _tbcMap;
  BatteryThermalVCMap _tvcMap;

  BatteryModelOptions<T> _options;

  MFRPtr _initialPotentialNorm;
  MFRPtr _initialThermalNorm;
  MFRPtr _initialPCNorm;
  vector<MFRPtr*> _initialSpeciesNormVector;
  int _niters;
  const int _nSpecies;

  BatteryModelFields _batteryModelFields;

  MFRPtr _currentPotentialResidual;
  MFRPtr _currentThermalResidual;
  MFRPtr _currentPCResidual;
  vector<MFRPtr*> _currentSpeciesResidual;
};

template<class T>
BatteryModel<T>::BatteryModel(GeomFields& geomFields,
			      const MeshList& realMeshes,
                              const MeshList& meshes,
                              const int nSpecies) :
  Model(meshes),
    _impl(new Impl(geomFields,realMeshes,meshes,nSpecies))
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
typename BatteryModel<T>::BatteryThermalBCMap&
BatteryModel<T>::getThermalBCMap() {return _impl->getThermalBCMap();}

template<class T>
typename BatteryModel<T>::BatteryThermalVCMap&
BatteryModel<T>::getThermalVCMap() {return _impl->getThermalVCMap();}

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
BatteryModel<T>::advanceThermal(const int niter)
{
  _impl->advanceThermal(niter);
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
void
BatteryModel<T>::recoverLastTimestep()
{
  _impl->recoverLastTimestep();
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
BatteryModel<T>::getHeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
  return _impl->getHeatFluxIntegral(mesh, faceGroupId);
}

template<class T>
T
BatteryModel<T>::getAverageConcentration(const Mesh& mesh, const int m)
{
  return _impl->getAverageConcentration(mesh, m);
}

template<class T>
T
BatteryModel<T>::getSpeciesResidual(const int speciesId)
{
  return _impl->getSpeciesResidual(speciesId);
}

template<class T>
T
BatteryModel<T>::getPotentialResidual()
{
  return _impl->getPotentialResidual();
}

template<class T>
T
BatteryModel<T>::getThermalResidual()
{
  return _impl->getThermalResidual();
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

template<class T>
T
BatteryModel<T>::getFaceGroupArea(const Mesh& mesh, const int fgID)
{
  return _impl->getFaceGroupArea(mesh, fgID);
}

template<class T>
T
BatteryModel<T>::getFaceGroupVoltage(const Mesh& mesh, const int fgID)
{
  return _impl->getFaceGroupVoltage(mesh, fgID);
}

template<class T>
T
BatteryModel<T>::getMeshVolume(const Mesh& mesh)
{
  return _impl->getMeshVolume(mesh);
}

