// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
#include "ElecDiffusionDiscretization.h"
#include "DriftDiscretization.h"
#include "TimeDerivativeDiscretization.h"
#include "TunnelingDiscretization.h"
#include "EmissionDiscretization.h"
#include "CaptureDiscretization.h"
#include "InjectionDiscretization.h"
#include "TrapBandTunnelingDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "SourceDiscretization.h"

#include "PhysicsConstant.h"
#include "ElectricUtilityFunctions.h"

#include "SquareTensor.h"
#include "ElecDiagonalTensor.h"
//#include "LinearizeDielectric.h"
#include "LinearizeInterfaceJump.h"
#include "LinearizePotentialInterface.h"
#include "BatteryElectricDiffusionDiscretization.h"

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

template<class T>
class ElectricModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Array<int> IntArray;
  typedef Vector<T,3> VectorT3;
  typedef Vector<double, 3> VectorD3;
  typedef Array<VectorT3> VectorT3Array;
  
  typedef Vector<T,3> VectorTN;
  typedef Array<VectorTN> VectorTNArray;
  
  typedef SquareTensor<T, 3> TensorNxN;
  //typedef ElecDiagonalTensor<T, 2> TensorNxN;
  typedef Array<TensorNxN> TensorNxNArray;

  //typedef ElecOffDiagonalTensor<T, 2> OffDiag;

  typedef Gradient<VectorTN> CGradType;
  typedef Array<Gradient<VectorTN> > CGradArray;    
  
  typedef Gradient<T> PGradType;
  typedef Array<Gradient<T> > PGradArray;
  
  

  Impl(const GeomFields& geomFields,
       ElectricFields& electricFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _electricFields(electricFields),
    _potentialGradientModel(_meshes,_electricFields.potential,
                              _electricFields.potential_gradient,_geomFields),
    _chargeGradientModel (_meshes, _electricFields.charge,
			  _electricFields.chargeGradient,_geomFields),
    _initialElectroStaticsNorm(),
    _initialChargeTransportNorm(),
    _niters(0),
    _avgCharge(0),
    _tunnelCurrentIn(0),
    _tunnelCurrentOut(0)
  
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

	ElectricVC<T> *vc(new ElectricVC<T>());
        vc->vcType = "dielectric";
        _vcMap[mesh.getID()] = vc;
        
 	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            ElectricBC<T> *bc(new ElectricBC<T>());
            
            _bcMap[fg.id] = bc;
            if (fg.groupType == "wall") 
            {
                bc->bcType = "SpecifiedPotential";
            }
	    else if (fg.groupType == "symmetry") 
            {
                bc->bcType = "Symmetry";
            }
	    //else
            //  throw CException("ElectricModel: unknown face group type "
            //                   + fg.groupType);
        }
    }
  }

  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
	
	const ElectricVC<T>& vc = *_vcMap[mesh.getID()];

	const StorageSite& cells = mesh.getCells();

	const StorageSite& faces = mesh.getFaces();
	
	// initialize fields for electrostatics

	if (_options.electrostatics_enable){

	  const int nCells = cells.getCountLevel1();
	
	  //initial potential setup
	  shared_ptr<TArray> pCell(new TArray(nCells));
	  *pCell = _options["initialPotential"];
	  _electricFields.potential.addArray(cells,pCell);

	  //dielectric_constant setup
	  shared_ptr<TArray> permCell(new TArray(nCells));
	  *permCell = vc["dielectric_constant"] * E0_SI ;
	  //*permCell = vc["dielectric_constant"];
	  _electricFields.dielectric_constant.addArray(cells,permCell);
	  
	  //total_charge (source in Poisson equation) setup
	  if ( vc.vcType == "dielectric"){
	    shared_ptr<TArray> sdCell(new TArray(nCells));
	    *sdCell = _options["initialTotalCharge"]*-QE;
	    _electricFields.total_charge.addArray(cells,sdCell);
	  }
	  else{
	    shared_ptr<TArray> saCell(new TArray(nCells));
	    saCell->zero();
	    _electricFields.total_charge.addArray(cells,saCell);
	  }
	  //potential gradient setup
	  shared_ptr<PGradArray> gradp(new PGradArray(nCells));
          gradp->zero();	
          _electricFields.potential_gradient.addArray(cells,gradp);

	  //electric_field setup
	  shared_ptr<VectorT3Array> ef(new VectorT3Array(nCells));
          ef->zero();	
          _electricFields.electric_field.addArray(cells,ef);

	  //initial potential flux; Note: potential_flux only stored on boundary faces
	  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<TArray> pFlux(new TArray(faces.getCount()));
	      pFlux->zero();
	      _electricFields.potential_flux.addArray(faces,pFlux);
	      
	    }
	  foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<TArray> pFlux(new TArray(faces.getCount()));
	      pFlux->zero();
	      _electricFields.potential_flux.addArray(faces,pFlux);
	      
	    }

	  //initial species concentration setup for linking to SpeciesModel
	  shared_ptr<TArray> sCCell(new TArray(nCells));
	  sCCell->zero();
	  _electricFields.speciesConcentration.addArray(cells,sCCell);

	  shared_ptr<TArray> lnSCCell(new TArray(nCells));
	  lnSCCell->zero();
	  _electricFields.lnSpeciesConcentration.addArray(cells,lnSCCell);

	  shared_ptr<PGradArray> gradLnSC(new PGradArray(nCells));
	  gradLnSC->zero();
	  _electricFields.lnSpeciesConcentrationGradient.addArray(cells,gradLnSC);

	}

	//initilize fields in charge transport models
	/* 
	   regarding multi trapdepth model, assume the number of trapdepth is nTrap
	   index i = 0 to i = nTrap-1 refer to charges in traps
	   index i = nTrap refer to charges in conduction band
	*/

	if (_options.chargetransport_enable){

	  if ( vc.vcType == "dielectric" ) {

	    const int nCells = cells.getCountLevel1();

	    //conduction_band setup
	    shared_ptr<TArray> cb(new TArray(nCells));
	    cb->zero();
	    _electricFields.conduction_band.addArray(cells,cb);

	    //valence_band setup
	    shared_ptr<TArray> vb(new TArray(nCells));
	    vb->zero();
	    _electricFields.valence_band.addArray(cells, vb);

	    //charge density setup 
	    shared_ptr<VectorTNArray> ch(new VectorTNArray (nCells));
	    ch->zero();
	    _electricFields.charge.addArray(cells, ch);	    

	    if (_options.transient_enable)
	    {
	      _electricFields.chargeN1.addArray(cells,dynamic_pointer_cast<ArrayBase>(ch->newCopy()));
	      if (_options.timeDiscretizationOrder > 1)
		_electricFields.chargeN2.addArray(cells, dynamic_pointer_cast<ArrayBase>(ch->newCopy()));
	    }

	    //charge velocity
	    shared_ptr<VectorT3Array> vcell(new VectorT3Array(nCells));
	    vcell->zero();
	    _electricFields.electron_velocity.addArray(cells, vcell);
	    
	    //free_electron_capture_cross setup
	    shared_ptr<VectorTNArray> fecr(new VectorTNArray(nCells));
	    fecr->zero();
	    _electricFields.free_electron_capture_cross.addArray(cells, fecr);

	    // convection flux on all faces 
	    shared_ptr<TArray> mf (new TArray(faces.getCount()));
	    mf->zero();
	    _electricFields.convectionFlux.addArray(faces, mf);

	    //diffusivity 
	    shared_ptr<TArray> diffCell(new TArray(cells.getCountLevel1()));
	    const T diffCoeff = _constants["electron_mobility"] * K_SI * _constants["OP_temperature"] / QE;
	    *diffCell = diffCoeff;
	    _electricFields.diffusivity.addArray(cells,diffCell);
	
	    //create a zero field
	    shared_ptr<TArray> zeroCell(new TArray(cells.getCountLevel1()));
	    *zeroCell = T(0.0);
	    _electricFields.zero.addArray(cells,zeroCell);
	    
	    //create a one field
	    shared_ptr<TArray> oneCell(new TArray(cells.getCountLevel1()));
	    *oneCell = T(1.0);
	    _electricFields.one.addArray(cells,oneCell);

	    //initial charge gradient array
	    shared_ptr<CGradArray> gradC(new CGradArray(cells.getCountLevel1()));
	    gradC->zero();
	    _electricFields.chargeGradient.addArray(cells,gradC);
        
	    //initial charge flux
	    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<VectorTNArray> cFlux(new VectorTNArray(faces.getCount()));
	      cFlux->zero();
	      _electricFields.chargeFlux.addArray(faces,cFlux);
	    }
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<VectorTNArray> cFlux(new VectorTNArray(faces.getCount()));
	      cFlux->zero();
	      _electricFields.chargeFlux.addArray(faces,cFlux);
	    }
	  }
	}
    }
    
    _electricFields.dielectric_constant.syncLocal();
    _niters  = 0;
    _initialElectroStaticsNorm = MFRPtr();
    if (_options.chargetransport_enable)
      _initialChargeTransportNorm = MFRPtr();
  }
  

 

  ElectricBCMap& getBCMap() {return _bcMap;}

  ElectricBC<T>& getBC(const int id) {return *_bcMap[id];}

  ElectricVCMap& getVCMap() {return _vcMap;}

  ElectricVC<T>& getVC(const int id) {return *_vcMap[id];}
  
  ElectricModelOptions<T>& getOptions() {return _options;}

  ElectricModelConstants<T>& getConstants() {return _constants;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)    {
      
        const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();
	
	VectorTNArray& charge =
          dynamic_cast<VectorTNArray&>(_electricFields.charge[cells]);
        VectorTNArray& chargeN1 =
          dynamic_cast<VectorTNArray&>(_electricFields.chargeN1[cells]);
	TArray& totalcharge = 
	  dynamic_cast<TArray&> (_electricFields.total_charge[cells]);

        if (_options.timeDiscretizationOrder > 1)
        {
            VectorTNArray& chargeN2 =
              dynamic_cast<VectorTNArray&>(_electricFields.chargeN2[cells]);
            chargeN2 = chargeN1;
        }
        chargeN1 = charge;
	/*
	const int nTrap = _constants["nTrap"];
	for (int c=0; c<nCells; c++){
	  for (int i=0; i<=nTrap; i++){
	    totalcharge[c] += charge[c][i]*-QE;
	  }	  
	}
	*/
	const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
	if (vc.vcType == "dielectric"){
	  const TArray& cellVolume = 
	    dynamic_cast<const TArray&> (_geomFields.volume[cells]);
	  T chargeSum(0);
	  T volumeSum(0);
	  const int nCells = cells.getSelfCount();
	  for (int c=0; c<nCells; c++){
	    volumeSum += cellVolume[c];
	    chargeSum += totalcharge[c]*cellVolume[c];
	  }
	  _avgCharge = chargeSum/volumeSum;
	}
	cout << "the avergae charge " << _avgCharge << endl;
     }

  }


  MFRPtr solveElectroStatics()
  {
    LinearSystem ls;
        
    initElectroStaticsLinearization(ls);
    
    ls.initAssembly();
   
    linearizeElectroStatics(ls);
   
    ls.initSolve();
   
    MFRPtr rNorm(_options.getElectroStaticsLinearSolver().solve(ls));
  
    if (!_initialElectroStaticsNorm) _initialElectroStaticsNorm = rNorm;
   
    _options.getElectroStaticsLinearSolver().cleanup();
   
    ls.postSolve();
  
    ls.updateSolution();
  
    updateElectricField();
  
    if (_options.chargetransport_enable){

      updateElectronVelocity();
      
      updateConvectionFlux();
    }
  
    return rNorm;
   
  }

  MFRPtr solveChargeTransport()
  {
    LinearSystem ls;
        
    initChargeTransportLinearization(ls);

    ls.initAssembly();

    linearizeChargeTransport(ls);

    ls.initSolve();

    MFRPtr rNorm(_options.getChargeTransportLinearSolver().solve(ls));

    if (!_initialChargeTransportNorm) _initialChargeTransportNorm = rNorm;
     
    _options.getChargeTransportLinearSolver().cleanup();

    ls.postSolve();

    ls.updateSolution();

    
    return rNorm;
  }

  void initElectroStaticsLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex tIndex(&_electricFields.potential,&cells);

        ls.getX().addArray(tIndex,_electricFields.potential.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(tIndex,tIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_electricFields.potential_flux,&faces);
            ls.getX().addArray(fIndex,_electricFields.potential_flux.getArrayPtr(faces));

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

            MultiField::ArrayIndex fIndex(&_electricFields.potential_flux,&faces);
            ls.getX().addArray(fIndex,_electricFields.potential_flux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,tIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

    }
  }

   void initChargeTransportLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

	//const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
	
	//if (vc.vcType == "dielectric"){
        
	  const StorageSite& cells = mesh.getCells();
	
	  MultiField::ArrayIndex cIndex(&_electricFields.charge,&cells);

	  ls.getX().addArray(cIndex,_electricFields.charge.getArrayPtr(cells));

	  const CRConnectivity& cellCells = mesh.getCellCells();
        
	  shared_ptr<Matrix> m(new CRMatrix<TensorNxN,TensorNxN,VectorTN>(cellCells));

	  ls.getMatrix().addMatrix(cIndex,cIndex,m);

	  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;

	      MultiField::ArrayIndex fIndex(&_electricFields.chargeFlux,&faces);
	      ls.getX().addArray(fIndex,_electricFields.chargeFlux.getArrayPtr(faces));

	      const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	      shared_ptr<Matrix> mft(new FluxJacobianMatrix<TensorNxN,VectorTN>(faceCells));
	      ls.getMatrix().addMatrix(fIndex,cIndex,mft);

	      shared_ptr<Matrix> mff(new DiagonalMatrix<TensorNxN,VectorTN>(faces.getCount()));
	      ls.getMatrix().addMatrix(fIndex,fIndex,mff);
	    }

	  foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      MultiField::ArrayIndex fIndex(&_electricFields.chargeFlux,&faces);
	      ls.getX().addArray(fIndex,_electricFields.chargeFlux.getArrayPtr(faces));
	      
	      const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	      shared_ptr<Matrix> mft(new FluxJacobianMatrix<TensorNxN,VectorTN>(faceCells));
	      ls.getMatrix().addMatrix(fIndex,cIndex,mft);

	      shared_ptr<Matrix> mff(new DiagonalMatrix<TensorNxN,VectorTN>(faces.getCount()));
	      ls.getMatrix().addMatrix(fIndex,fIndex,mff);
	    }
	  //}
    }
  }

  void linearizeElectroStatics(LinearSystem& ls)
  {
    _potentialGradientModel.compute();
    
    DiscrList discretizations;
    
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>(_meshes,_geomFields,
                                            _electricFields.potential,
                                            _electricFields.dielectric_constant,
                                            _electricFields.potential_gradient,
					    _constants["dielectric_thickness"]));
    discretizations.push_back(dd);    

    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  _electricFields.potential,
	  _electricFields.total_charge
	));
    discretizations.push_back(sd);

    if(_options.ibm_enable){
      shared_ptr<Discretization>
	ibm(new GenericIBDiscretization<T,T,T>
	    (_meshes,_geomFields,_electricFields.potential));
      discretizations.push_back(ibm);
    }

    if(_options.ButlerVolmer)
      {
	//populate lnSpeciesConc Field
	for (int n=0; n<_meshes.size(); n++)
	  {
	    const Mesh& mesh = *_meshes[n];
	    const StorageSite& cells = mesh.getCells();
	    const TArray& speciesConc = dynamic_cast<const TArray&>(_electricFields.speciesConcentration[cells]);
	    TArray& lnSpeciesConc = dynamic_cast<TArray&>(_electricFields.lnSpeciesConcentration[cells]);
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
	GradientModel<T> lnSpeciesGradientModel(_meshes,_electricFields.lnSpeciesConcentration,
					  _electricFields.lnSpeciesConcentrationGradient,_geomFields);
	lnSpeciesGradientModel.compute();
	
	//discretize ln term (only affects electrolyte mesh)
	shared_ptr<Discretization>
	  bedd(new BatteryElectricDiffusionDiscretization<T,T,T>(_meshes,_geomFields,
						_electricFields.potential,
						_electricFields.dielectric_constant,
						_electricFields.lnSpeciesConcentration,
						_electricFields.lnSpeciesConcentrationGradient,
						_options["BatteryElectrolyteMeshID"]));
	discretizations.push_back(bedd);    

      }
    
    Linearizer linearizer;

    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
    
    

    
    const int numMeshes = _meshes.size();

    /*linearize shell mesh*/
    /*
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      if (mesh.isShell())
	{	    
	  LinearizeDielectric<T, T, T> lsm (_geomFields,
					    _electricFields.dielectric_constant,
					    _constants["dielectric_thickness"],
					    _electricFields.potential);

	  lsm.discretize(mesh, ls.getMatrix(), ls.getX(), ls.getB() );
	}
    }
    */

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
	      LinearizePotentialInterface<T, T, T> lbv (_geomFields,
							_electricFields.potential,
							_electricFields.speciesConcentration,
							_options["ButlerVolmerRRConstant"],
							_options["Interface_A_coeff"],
							_options["Interface_B_coeff"],
							Anode,
							Cathode);

	      lbv.discretize(mesh, parentMesh, otherMesh, ls.getMatrix(), ls.getX(), ls.getB() );
	    }
	  else
	    {
	  
	      LinearizeInterfaceJump<T, T, T> lsm (_options["Interface_A_coeff"],
						   _options["Interface_B_coeff"],
						   _electricFields.potential);

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
            const ElectricBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _electricFields.potential,
                                  _electricFields.potential_flux,
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
	    else if (bc.bcType == "SpecialDielectricBoundary")
	    {
	    	FloatValEvaluator<T> bT(bc.getVal("specifiedPotential"), faces);
	    	//const T specifiedPotential(bc["specifiedPotential"]);
	    	const T dielectric_thickness = _constants["dielectric_thickness"];
	    	const T dielectric_constant = _constants["dielectric_constant"] * E0_SI;
	    	const T coeff = dielectric_constant / dielectric_thickness;
	    	//const T avgCharge = 100/1.60217646e-19*8.854187826e-12*-QE;
	    	_avgCharge = -1e4;
	    	//const T src = _avgCharge * dielectric_thickness;
	    	const T src = 0;
	    	for(int f=0; f<nFaces; f++)
			gbc.applyDielectricInterfaceBC(f, coeff, bT[f], src);
	    }
	    else
                throw CException(bc.bcType + " not implemented for ElectricModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
	  
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

	    
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _electricFields.potential,
                                  _electricFields.potential_flux, 
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }
  


 void linearizeChargeTransport(LinearSystem& ls)
  {
    _chargeGradientModel.compute();
    
    DiscrList discretizations;
    
    if (_options.tunneling_enable){
      shared_ptr<Discretization>
	tnd(new TunnelingDiscretization<VectorTN, TensorNxN, TensorNxN>
	    (_meshes, _geomFields,
	     _electricFields.charge,	     
	     _electricFields.conduction_band,
	     _constants, 
	     _tunnelCurrentIn, 
	     _tunnelCurrentOut
	     ));
      discretizations.push_back(tnd);      
    }
    
    if (_options.injection_enable){
      shared_ptr<Discretization>
	inj(new InjectionDiscretization<VectorTN, TensorNxN, TensorNxN>
	    (_meshes, _geomFields,
	     _electricFields.charge,
	     _electricFields.electric_field,
	     _electricFields.conduction_band,
	     _constants
	     ));
      discretizations.push_back(inj);
    }
    
    if (_options.emission_enable){
      shared_ptr<Discretization>
	em(new EmissionDiscretization<VectorTN, TensorNxN, TensorNxN>
	   (_meshes, _geomFields,
	    _electricFields.charge,
	    _electricFields.electric_field,
	    _constants));
      discretizations.push_back(em);
      }
    
    if (_options.capture_enable){
      shared_ptr<Discretization>
	capt(new CaptureDiscretization<VectorTN, TensorNxN, TensorNxN>
	   (_meshes, _geomFields,
	    _electricFields.charge,
	    _electricFields.free_electron_capture_cross,
	    _constants));
      discretizations.push_back(capt);
    }
   
    
    if (_options.trapbandtunneling_enable){
      shared_ptr<Discretization>
	tbt(new TrapBandTunnelingDiscretization<VectorTN, TensorNxN, TensorNxN>
	    (_meshes, _geomFields,
	     _electricFields.charge,
	     _electricFields.electric_field,
	     _electricFields.conduction_band,
	     _constants));
      discretizations.push_back(tbt);

    }
    
#if 0
    if (_options.diffusion_enable){
      shared_ptr<Discretization>
	dd(new ElecDiffusionDiscretization<VectorTN, TensorNxN, TensorNxN>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.diffusivity,
	    _electricFields.chargeGradient));
      discretizations.push_back(dd);
    }
#endif

    if (_options.drift_enable){
      shared_ptr<Discretization>
	cd(new DriftDiscretization<VectorTN, TensorNxN, TensorNxN>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.convectionFlux,
	    _constants["nTrap"]));
      discretizations.push_back(cd);
    }
    
    if (_options.transient_enable)
      {
      shared_ptr<Discretization>
	td(new TimeDerivativeDiscretization<VectorTN, TensorNxN, TensorNxN>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.chargeN1,
	    _electricFields.chargeN2,
	    _electricFields.one,
	    _options["timeStep"]));
      discretizations.push_back(td);
      }

   

    
    Linearizer linearizer;
      
    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());
   
    if (_options.diffusion_enable || _options.drift_enable){
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            const ElectricBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<VectorTN,TensorNxN,TensorNxN> gbc(faces,mesh,
                                  _geomFields,
                                  _electricFields.charge,
                                  _electricFields.chargeFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());
	    //dielectric charging uses fixed zero dirichlet bc
	    VectorTN zero(VectorTN::getZero());;
	        
	    //gbc.applyNonzeroDiagBC();
	    if (bc.bcType == "Symmetry")
	      //gbc.applyExtrapolationBC();
	      gbc.applyDirichletBC(zero);
	    else
	      gbc.applyDirichletBC(zero);
	    
	  }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<VectorTN,TensorNxN,TensorNxN> gbc(faces,mesh,
                                  _geomFields,
                                  _electricFields.charge,
                                  _electricFields.chargeFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());
            gbc.applyInterfaceBC();
	  }
	}
    }
   
    
  }



 
  bool advance(const int niter)
  {
    bool flag1 = false;
    bool flag2 = false;
         
    if(_options.electrostatics_enable){
     
      for(int n=0; n<niter; n++)
	{
	  MFRPtr eNorm = solveElectroStatics();
	  	
	  if (_niters < 5)
	    _initialElectroStaticsNorm->setMax(*eNorm);
	  
	  MFRPtr eNormRatio(eNorm->normalize(*_initialElectroStaticsNorm));

#ifndef FVM_PARALLEL
	  if (_options.printNormalizedResiduals)
	    cout << _niters << ": " << *eNormRatio << ";" <<  endl;
	  else
	    cout << _niters << ": " << *eNorm << ";" <<  endl;
#endif

#ifdef FVM_PARALLEL
     if ( MPI::COMM_WORLD.Get_rank() == 0 ){
	  if (_options.printNormalizedResiduals)
	    cout << _niters << ": " << *eNormRatio << ";" <<  endl;
	  else
	    cout << _niters << ": " << *eNorm << ";" <<  endl;
     }
#endif
   
	  if (*eNormRatio < _options.electrostaticsTolerance ) {
	      flag1 = true;
	      break;
	  }

	  
	}
    }

      
    if(_options.chargetransport_enable){
      for(int n=0; n<niter; n++)
	{
	  generateBandDiagram();
	  
	  MFRPtr cNorm = solveChargeTransport();
	  
	  if (_niters < 5)
	    _initialChargeTransportNorm->setMax(*cNorm);
	  
	  MFRPtr cNormRatio((*cNorm)/(*_initialChargeTransportNorm));
        
	  if (_options.printNormalizedResiduals)
	    cout << _niters << ": " << *cNormRatio <<  endl;
	  else
	    cout << _niters << ": " <<  *cNorm <<  endl;

	  if (*cNormRatio < _options.chargetransportTolerance){
	    flag2 = true;
	    //break;
	  }
	}

    }	
    _niters++;  
    return flag1;
   
  }

  /*** update electric_field from potential_gradient ***/
  void updateElectricField()
  {
    _potentialGradientModel.compute();

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      VectorT3Array& electric_field = dynamic_cast<VectorT3Array& > (_electricFields.electric_field[cells]);
      const PGradArray& potential_gradient = dynamic_cast<const PGradArray& > (_electricFields.potential_gradient[cells]);

      for(int c=0; c<nCells; c++){
	electric_field[c][0] = -potential_gradient[c][0];
	electric_field[c][1] = -potential_gradient[c][1];
	electric_field[c][2] = -potential_gradient[c][2];
      }
    }
  }

  /*** update electron velocity after electric_field is updated ***/
  void updateElectronVelocity()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCountLevel1();
      const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array& > (_electricFields.electric_field[cells]);
      VectorT3Array& electron_velocity = dynamic_cast<VectorT3Array& > (_electricFields.electron_velocity[cells]);
      const T electron_mobility = _constants["electron_mobility"];
      const T electron_saturation_velocity =  _constants["electron_saturation_velocity"];

      // need to convert to 3D 
      for(int c=0; c<nCells; c++){
	VectorT3 vel = electron_mobility * electric_field[c];
	if (mag(vel) < electron_saturation_velocity ){
	  electron_velocity[c] = - electron_mobility * electric_field[c];
	}
	else {
	  electron_velocity[c] = -  electron_saturation_velocity * (electric_field[c] / mag(electric_field[c]));
	}
      }     
    }
   
  }

  void updateConvectionFlux()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      const StorageSite& faces = mesh.getFaces();
      const int nFaces = faces.getCount();
      const VectorT3Array& vel = dynamic_cast<const VectorT3Array& > (_electricFields.electron_velocity[cells]);
      TArray& convFlux = dynamic_cast<TArray&> (_electricFields.convectionFlux[faces]);
      const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
      const CRConnectivity& faceCells = mesh.getAllFaceCells();
      
      for ( int f=0; f<nFaces; f++){
	const int c0 = faceCells(f, 0);
	const int c1 = faceCells(f, 1);
	convFlux[f] = 0.5 * (dot(vel[c0], faceArea[f]) + dot(vel[c1], faceArea[f]));	
      }
      
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  TArray& convFlux = dynamic_cast<TArray&> (_electricFields.convectionFlux[faces]);
	  const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	  const CRConnectivity& faceCells = mesh.getAllFaceCells();
	  const ElectricBC<T>& bc = *_bcMap[fg.id];
	  if(bc.bcType == "Symmetry")
	    {
	      convFlux.zero();
	    }
	  else{
	    for(int f=0; f<nFaces; f++){
	      const int c0 = faceCells(f,0);
	      convFlux[f] = dot(vel[c0], faceArea[f]);	     
	    }
	    
	  }
	}      
    }
  }
      
  /*** Generate Band Diagram after Poisson solve is done ***/
  void generateBandDiagram()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
	
      const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
        
      const StorageSite& cells = mesh.getCells(); 
      
      const int nCells = cells.getCountLevel1();
      
      const T& dielectric_ionization = _constants["dielectric_ionization"];
      const T& dielectric_bandgap = _constants["dielectric_bandgap"];
      const TArray& potential = dynamic_cast<const TArray& > (_electricFields.potential[cells]);
      
      TArray& conduction_band = dynamic_cast< TArray& > (_electricFields.conduction_band[cells]);
      TArray& valence_band = dynamic_cast< TArray& > (_electricFields.valence_band[cells]);
      // ??? what is 2 and 4 //
      if (vc.vcType == "dielectric"){
	for(int c=0; c<nCells; c++){
	  conduction_band[c] = -(dielectric_ionization + potential[c]);
	  valence_band[c] = conduction_band[c] - dielectric_bandgap;
	}
      }
      else if (vc.vcType == "air"){
	for(int c=0; c<nCells; c++){
	  conduction_band[c] = -(dielectric_ionization + potential[c]) + 2;
	  valence_band[c] = conduction_band[c] - dielectric_bandgap - 4;
	}
      }
      else 
	throw CException(vc.vcType + " not implemented for ElectricModel");
    }
  }


  void calculateEquilibriumParameters()
  {

    /* this gives the real initial condition for charges */

    const T membrane_workfunction = _constants["membrane_workfunction"];
    const T substrate_workfunction = _constants["substrate_workfunction"];
    const T dielectric_thickness = _constants["dielectric_thickness"];
    const T optical_dielectric_constant = _constants["optical_dielectric_constant"];
    const T dielectric_ionization = _constants["dielectric_ionization"];   
    const T temperature = _constants["OP_temperature"];
    const T electron_effmass = _constants["electron_effmass"];
    const T poole_frenkel_emission_frequency = _constants["poole_frenkel_emission_frequency"];
    const int normal = _constants["normal_direction"];
    const int nTrap = _constants["nTrap"];
    vector<T> electron_trapdensity = _constants.electron_trapdensity;
    vector<T> electron_trapdepth = _constants.electron_trapdepth;

    if (int(electron_trapdensity.size()) != nTrap || int(electron_trapdepth.size()) != nTrap){
      throw CException ("wrong trapdepth size!");
    }

    T effefield = (membrane_workfunction - substrate_workfunction) / dielectric_thickness;

    T alpha = sqrt(QE / (PI * E0_SI * optical_dielectric_constant));

    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
      if (vc.vcType == "dielectric"){

	const StorageSite& cells = mesh.getCells(); 
	const int nCells = cells.getCountLevel1();

	const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);

	VectorTNArray& free_electron_capture_cross = dynamic_cast<VectorTNArray&> (_electricFields.free_electron_capture_cross[cells]);	

	VectorTNArray& charge = dynamic_cast<VectorTNArray& > (_electricFields.charge[cells]);

 	VectorTNArray& chargeN1 = dynamic_cast<VectorTNArray& > (_electricFields.chargeN1[cells]);

        VectorTNArray* chargeN2  = 0;

	if (_options.timeDiscretizationOrder > 1){
	  chargeN2 = &(dynamic_cast<VectorTNArray& > (_electricFields.chargeN2[cells]));
	}

	for(int c=0; c<nCells; c++)
	{ 	  
	  
	  T fermilevel = -substrate_workfunction+effefield*cellCentroid[c][normal];  
	  T energy(0.);	  
	  
	  for (int i=0; i<nTrap; i++){
	    
	    energy = -dielectric_ionization - electron_trapdepth[i];

	    charge[c][i] = electron_trapdensity[i] * FermiFunction(energy, fermilevel, temperature);
	    chargeN1[c][i] = charge[c][i];
	  
	    if (_options.timeDiscretizationOrder > 1)
	      (*chargeN2)[c][i] = charge[c][i];
	  
	    energy = -dielectric_ionization;
	  
	    charge[c][nTrap] += electron_trapdensity[i] * FermiFunction(energy, fermilevel, temperature);
	    chargeN1[c][nTrap] = charge[c][nTrap];
	  

	    if (_options.timeDiscretizationOrder > 1)
	      (*chargeN2)[c][nTrap] = charge[c][nTrap];
	  }

	  for (int i=0; i<nTrap; i++){
	    T expt = (electron_trapdepth[i] - alpha * sqrt(fabs(effefield))) * QE / (K_SI*temperature);
	    T beta = exp(-expt);
	    T velocity = sqrt(8 * K_SI * temperature / (PI * ME * electron_effmass));
	    
	    free_electron_capture_cross[c][i] = charge[c][i] * poole_frenkel_emission_frequency * beta;
	    free_electron_capture_cross[c][i] /= (velocity*(electron_trapdensity[i]-charge[c][i])*charge[c][nTrap]);
	    
	  }
	}
      }
    }
    
  }

 
  void computeIBFacePotential(const StorageSite& solid)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    
    const TArray& sP =
      dynamic_cast<const TArray&>(_electricFields.potential[solid]);

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
	if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0)
        {
        	const StorageSite& cells = mesh.getCells();
        	const StorageSite& ibFaces = mesh.getIBFaces();
        
        	GeomFields::SSPair key1(&ibFaces,&cells);
        	const IMatrix& mIC =
          	dynamic_cast<const IMatrix&>
          	(*_geomFields._interpolationMatrices[key1]);

		GeomFields::SSPair key2(&ibFaces,&solid);
        	const IMatrix& mIP =
          	dynamic_cast<const IMatrix&>
          	(*_geomFields._interpolationMatrices[key2]);

		shared_ptr<TArray> ibP(new TArray(ibFaces.getCount()));
        
        	const TArray& cP =
          	dynamic_cast<const TArray&>(_electricFields.potential[cells]);

		//const Array<T>& cellToIBCoeff = mIC.getCoeff();
		//const Array<T>& particlesToIBCoeff = mIP.getCoeff();
	
		ibP->zero();

		mIC.multiplyAndAdd(*ibP,cP);
		mIP.multiplyAndAdd(*ibP,sP);

        	_electricFields.potential.addArray(ibFaces,ibP);
        	
        	
        }
    }

  }      


  void
  computeSolidSurfaceForce(const StorageSite& solidFaces, bool perUnitArea)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;

    const int nSolidFaces = solidFaces.getCount();

    updateElectricField();

    boost::shared_ptr<VectorT3Array>
      forcePtr( new VectorT3Array(nSolidFaces));
    VectorT3Array& force = *forcePtr;

    force.zero();
    _electricFields.force.addArray(solidFaces,forcePtr);

    const VectorT3Array& solidFaceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[solidFaces]);

    const TArray& solidFaceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[solidFaces]);

    //const VectorT3Array& solidFaceCoordinate =
    //  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[solidFaces]);

       
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        if (mesh.isShell())
        	return;
        const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
        const StorageSite& cells = mesh.getCells();

      	const VectorT3Array& electric_field = 
		dynamic_cast<const VectorT3Array&> (_electricFields.electric_field[cells]);

      	const TArray& dielectric_constant = 
		dynamic_cast<const TArray&>(_electricFields.dielectric_constant[cells]);
      
      	const CRConnectivity& solidFacesToCells
		= mesh.getConnectivity(solidFaces,cells);
        
      	const IntArray& sFCRow = solidFacesToCells.getRow();
      	const IntArray& sFCCol = solidFacesToCells.getCol();
      
      	GeomFields::SSPair key1(&solidFaces,&cells);
      	const IMatrix& mIC = dynamic_cast<const IMatrix&>
		(*_geomFields._interpolationMatrices[key1]);

      	const Array<T>& iCoeffs = mIC.getCoeff();

      	for(int f=0; f<solidFaces.getCount(); f++)
      	{
	  T forceMag(0);
	  T forceSign(1);
	  //cout << f << endl;
	  for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)	    {
	    const int c = sFCCol[nc];
	    const T coeff = iCoeffs[nc];
	    T efmag = mag(electric_field[c]);
	    
	    forceSign = dot(electric_field[c], solidFaceArea[f]);
	    
	    if (fabs(forceSign) > 0.0) 
	      forceSign /= fabs(forceSign);
	    else{
	      forceSign = 0.0;
	    }
	    //T coeff = 0.25;
	    forceMag += 0.5 * coeff* dielectric_constant[c] *  efmag * efmag * forceSign; 
	    //cout << "   " << c << " " << efmag << " "  << coeff << " " << forceSign << endl;
	  }
	  //cout << forceMag << endl;

	  const VectorT3& Af = solidFaceArea[f];
	  force[f] += Af * forceMag;
	 
	  if (perUnitArea){
	    force[f] /= solidFaceAreaMag[f];
	  }

	}
      	  
    }
  }

  void printBCs()
  {
    foreach(typename ElectricBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename ElectricBC<T>::value_type& vp, *pos.second)
        {
            cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
        }
    }
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
              dynamic_cast<const TArray&>(_electricFields.potential_flux[faces]);
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
              dynamic_cast<const TArray&>(_electricFields.potential_flux[faces]);
            for(int f=0; f<nFaces; f++)
              r += potentialFlux[f];
            found=true;
	  }
    }


    if (!found)
      throw CException("getPotentialFluxIntegral: invalid faceGroupID");
    return r;
  }

  map<string,shared_ptr<ArrayBase> >&
  getPersistenceData()
  {
    _persistenceData.clear();
    
    Array<int>* niterArray = new Array<int>(1);
    (*niterArray)[0] = _niters;
    _persistenceData["niters"]=shared_ptr<ArrayBase>(niterArray);

    if (_initialElectroStaticsNorm)
    {
      _persistenceData["initialElectroStaticsNorm"] =
	_initialElectroStaticsNorm->getArrayPtr(_electricFields.potential);
    }
    else
    {
      Array<T>* xArray = new Array<T>(1);
      xArray->zero();
      _persistenceData["initialElectroStaticsNorm"]=shared_ptr<ArrayBase>(xArray);
    }
    return _persistenceData;
  }

  void restart()
  {
    if (_persistenceData.find("niters") != _persistenceData.end())
    {
        shared_ptr<ArrayBase> rp = _persistenceData["niters"];
        ArrayBase& r = *rp;
        Array<int>& niterArray = dynamic_cast<Array<int>& >(r);
        _niters = niterArray[0];
    }
    if (_persistenceData.find("initialElectroStaticsNorm") != _persistenceData.end())
    {
        shared_ptr<ArrayBase>  r = _persistenceData["initialElectroStaticsNorm"];
        _initialElectroStaticsNorm = MFRPtr(new MultiFieldReduction());
        _initialElectroStaticsNorm->addArray(_electricFields.potential,r);
    }
  }

  
  vector<T> getTunnelCurrent(){
    vector<T> flux;    
    flux.push_back(_tunnelCurrentIn);
    flux.push_back(_tunnelCurrentOut);
    return flux;
  }


    


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  ElectricFields& _electricFields;

  ElectricBCMap _bcMap;
  ElectricVCMap _vcMap;

  ElectricModelOptions<T> _options; 
  ElectricModelConstants<T> _constants;
  GradientModel<T> _potentialGradientModel;
  GradientModel<VectorTN> _chargeGradientModel;
  
  MFRPtr _initialElectroStaticsNorm;
  MFRPtr _initialChargeTransportNorm;
  int _niters;
  T _avgCharge;
  T _tunnelCurrentIn;
  T _tunnelCurrentOut;

 

  map<string,shared_ptr<ArrayBase> > _persistenceData;
};



template<class T>
ElectricModel<T>::ElectricModel(const GeomFields& geomFields,
                              ElectricFields& electricFields,
                              const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,electricFields,meshes))
{
  logCtor();
}


template<class T>
ElectricModel<T>::~ElectricModel()
{
  logDtor();
}

template<class T>
void
ElectricModel<T>::init()
{
  _impl->init();
}
 
/*
template<class T>
void 
ElectricModel<T>::dielectricOneDimModelPrep(const int nXCol, const int nYCol,const int nGrid,
				  const ElectricModel<T>::VectorD3 corner1_1, 
				  const ElectricModel<T>::VectorD3 corner1_2, 
				  const ElectricModel<T>::VectorD3 corner1_3,
				  const ElectricModel<T>::VectorD3 corner1_4,
				  const ElectricModel<T>::VectorD3 corner2_1, 
				  const ElectricModel<T>::VectorD3 corner2_2, 
				  const ElectricModel<T>::VectorD3 corner2_3,
				  const ElectricModel<T>::VectorD3 corner2_4)
{
  _impl->dielectricOneDimModelPrep (nXCol,  nYCol, nGrid,
				    corner1_1, corner1_2, corner1_3, corner1_4,
				    corner2_1, corner2_2, corner2_3, corner2_4) ;
}
*/

template<class T>
void
ElectricModel<T>::updateTime()
{
  _impl->updateTime();
}

template<class T>
typename ElectricModel<T>::ElectricBCMap&
ElectricModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
ElectricBC<T>&
ElectricModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
typename ElectricModel<T>::ElectricVCMap&
ElectricModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
ElectricVC<T>&
ElectricModel<T>::getVC(const int id) {return _impl->getVC(id);}

template<class T>
ElectricModelOptions<T>&
ElectricModel<T>::getOptions() {return _impl->getOptions();}

template<class T>
ElectricModelConstants<T>&
ElectricModel<T>::getConstants() {return _impl->getConstants();}

template<class T>
void
ElectricModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
bool
ElectricModel<T>::advance(const int niter)
{
  return _impl->advance(niter);  
}

template<class T>
void
ElectricModel<T>::calculateEquilibriumParameters()
{
  _impl->calculateEquilibriumParameters();
}

template<class T>
void
ElectricModel<T>::computeIBFacePotential(const StorageSite& particles)
{
  _impl->computeIBFacePotential(particles);
}
 
template<class T>
void
ElectricModel<T>::computeSolidSurfaceForce(const StorageSite& particles)
{
  return _impl->computeSolidSurfaceForce(particles,false);
}

template<class T>
T
ElectricModel<T>::getPotentialFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
  return _impl->getPotentialFluxIntegral(mesh, faceGroupId);
}

template<class T>
void
ElectricModel<T>::computeSolidSurfaceForcePerUnitArea(const StorageSite& particles)
{
  return _impl->computeSolidSurfaceForce(particles,true);
}

template<class T>
map<string,shared_ptr<ArrayBase> >&
ElectricModel<T>::getPersistenceData()
{
  return _impl->getPersistenceData();
}

template<class T>
void
ElectricModel<T>::restart()
{
  _impl->restart();
}


template<class T>
vector<T>
ElectricModel<T>::getTunnelCurrent()
{
  return _impl->getTunnelCurrent();
}
