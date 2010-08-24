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
#include "Octree.h"
#include "PhysicsConstant.h"
#include "ElectricUtilityFunctions.h"
#include "DielectricOneDimColumn.h"
#include "SquareTensor.h"


template<class T>
class ElectricModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Array<int> IntArray;
  typedef Vector<T,3> VectorT3;
  typedef Vector<double, 3> VectorD3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Vector<T,2> VectorT2;
  typedef Array<VectorT2> VectorT2Array;
  
  typedef SquareTensor<T, 2> Tensor2x2;
  typedef Array<Tensor2x2> Tensor2x2Array;
  
  typedef Gradient<T> PGradType;
  typedef Array<Gradient<T> > PGradArray;

  typedef Gradient<VectorT2> CGradType;
  typedef Array<Gradient<VectorT2> > CGradArray;
  
  typedef DielectricOneDimColumn<T> Column;
  typedef vector<shared_ptr<Column> > ColumnList;

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
    _columns(0), 
    _columnList(),
    _cellColumns(),
    _columnCells()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        
	//each mesh is initialized as type "dielectric"
	//need to reassigned with proper types in python 

	ElectricVC<T> *vc(new ElectricVC<T>());
        vc->vcType = "dielectric";
        _vcMap[mesh.getID()] = vc;
        
	//boundary condiction is initialized
	//need to be reassigned properly in python

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
	    else
              throw CException("ElectricModel: unknown face group type "
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
	
	const ElectricVC<T>& vc = *_vcMap[mesh.getID()];

	const StorageSite& cells = mesh.getCells();

	const StorageSite& faces = mesh.getFaces();
	
	if (_options.electrostatics_enable){

	  const int nCells = cells.getCount();
	
	  //initial potential setup
	  shared_ptr<TArray> pCell(new TArray(nCells));
	  *pCell = _options["initialPotential"];
	  _electricFields.potential.addArray(cells,pCell);

	  //dielectric_constant setup
	  shared_ptr<TArray> permCell(new TArray(nCells));
	  *permCell = vc["dielectric_constant"] * E0_SI ;
	  //*permCell = T(1.0);
	  _electricFields.dielectric_constant.addArray(cells,permCell);
	  
	  //total_charge (source in Poisson equation) setup
	  // there is no charge in air gap
	  if ( vc.vcType == "air"){
	    shared_ptr<TArray> saCell(new TArray(nCells));
	    saCell->zero();
	    _electricFields.total_charge.addArray(cells,saCell);
	  }
	  // there is a initial charge in dielectric
	  if ( vc.vcType == "dielectric"){
	    shared_ptr<TArray> sdCell(new TArray(nCells));
	    *sdCell = _options["initialTotalCharge"];
	    _electricFields.total_charge.addArray(cells,sdCell);
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
	}

	if (_options.chargetransport_enable){

	  /* the following fields only used in charge transport in dielectric*/
	  // if ( vc.vcType == "dielectric" ) {
	    const int nCells = cells.getCount();
	    //conduction_band setup
	    shared_ptr<TArray> cb(new TArray(nCells));
	    cb->zero();
	    _electricFields.conduction_band.addArray(cells,cb);

	    //valence_band setup
	    shared_ptr<TArray> vb(new TArray(nCells));
	    vb->zero();
	    _electricFields.valence_band.addArray(cells, vb);

	    //electron_totaltraps setup
	    shared_ptr<TArray> tt(new TArray(nCells));
	    *tt = _constants["electron_trapdensity"];
	    _electricFields.electron_totaltraps.addArray(cells, tt);

	    //charge density setup 
	    //VectorT2Array containing the electron density in trap and conduction band
	    shared_ptr<VectorT2Array> ch(new VectorT2Array(nCells));
	    ch->zero();
	    _electricFields.charge.addArray(cells, ch);
	    
	    //charge velocity
	    shared_ptr<VectorT3Array> vcell(new VectorT3Array(nCells));
	    vcell->zero();
	    _electricFields.electron_velocity.addArray(cells, vcell);

	    if (_options.transient_enable)
	    {
	      _electricFields.chargeN1.addArray(cells,
						  dynamic_pointer_cast<ArrayBase>(ch->newCopy()));
	      if (_options.timeDiscretizationOrder > 1)
		_electricFields.chargeN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(ch->newCopy()));
	    }
	   
	    //free_electron_capture_cross setup
	    shared_ptr<TArray> fecr(new TArray(nCells));
	    fecr->zero();
	    _electricFields.free_electron_capture_cross.addArray(cells, fecr);

	    
	    // convection flux on all faces 
	    shared_ptr<TArray> mf (new TArray(faces.getCount()));
	    mf->zero();
	    _electricFields.convectionFlux.addArray(faces, mf);

	    //diffusivity 
	    shared_ptr<TArray> diffCell(new TArray(cells.getCount()));
	    const T diffCoeff = _constants["electron_mobility"] * K_SI * _constants["OP_temperature"] / QE;
	    *diffCell = diffCoeff;
	    _electricFields.diffusivity.addArray(cells,diffCell);
	
	    //create a zero field
	    shared_ptr<TArray> zeroCell(new TArray(cells.getCount()));
	    *zeroCell = T(0.0);
	    _electricFields.zero.addArray(cells,zeroCell);
	    
	    //create a one field
	    shared_ptr<TArray> oneCell(new TArray(cells.getCount()));
	    *oneCell = T(1.0);
	    _electricFields.one.addArray(cells,oneCell);

	    //initial charge gradient array
	    shared_ptr<CGradArray> gradC(new CGradArray(cells.getCount()));
	    gradC->zero();
	    _electricFields.chargeGradient.addArray(cells,gradC);
        
	    //initial charge flux
	    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<VectorT2Array> cFlux(new VectorT2Array(faces.getCount()));
	      cFlux->zero();
	      _electricFields.chargeFlux.addArray(faces,cFlux);
	    }
	    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<VectorT2Array> cFlux(new VectorT2Array(faces.getCount()));
	      cFlux->zero();
	      _electricFields.chargeFlux.addArray(faces,cFlux);
	    }
	}
    }
    
    _electricFields.dielectric_constant.syncLocal();
    _niters  = 0;
    _initialElectroStaticsNorm = MFRPtr();
    if (_options.chargetransport_enable)
      _initialChargeTransportNorm = MFRPtr();
  }
  

  /*** set up the dielectric one-D model ***/
  /*** only applied to dielectric zone ***/

  void dielectricOneDimModelPrep
    (const int nXCol, const int nYCol, const int nGrid,
     const VectorD3 corner1_1, 
     const VectorD3 corner1_2,
     const VectorD3 corner1_3,
     const VectorD3 corner1_4,
     const VectorD3 corner2_1, 
     const VectorD3 corner2_2, 
     const VectorD3 corner2_3, 
     const VectorD3 corner2_4)
  {

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
      if (vc.vcType == "dielectric"){
	const StorageSite& cells = mesh.getCells();	
	const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[mesh.getCells()]);
	const int nCol = nXCol * nYCol;
	
	_columns.setCount(nCol);

	for(int n=0; n<nCol; n++){
	  shared_ptr<Column> cl(new Column(cells, nGrid));
	  _columnList.push_back(cl);
	}
	
	//set centerline start and end points and its normalized direction
	setCenterLines(nXCol, nYCol, 
		       corner1_1, corner1_2, corner1_3, corner1_4,
		       corner2_1, corner2_2, corner2_3, corner2_4, 
		       _columnList, _columns);   

	//setup connectivity between cells and columns
	setConnectivityCellColumns(cells, _columns, cellCentroid, _columnList,_cellColumns);	  
	_columnCells = _cellColumns->getTranspose(); 
	outputConnectivityCellColumns(mesh, cellCentroid, _cellColumns);

	//fill up the local celllist for each column
	setCellList(_columnList, _columnCells);

	//discretize center line for each column    
	discretizeCenterLine(_columnList);

	//setup connectivity between cell and grids for each column
	setConnectivityCellGrids(_columnList, cellCentroid);

	//calculate interpolation matrix for each column
	const int method = 2;
	calculateInterpolationMatrix(_columnList, method, cellCentroid);

	/*	
	//test interpolation
	TArray& transmission = dynamic_cast<TArray&> (_electricFields.transmission[cells]);
	TArray& conduction_band = dynamic_cast< TArray&> (_electricFields.conduction_band[cells]);
	
	//test 1: assume conduction_band is linear with coordinate z;
	for(int c=0; c<cells.getCount(); c++){
	  conduction_band[c] = (cellCentroid[c][2]) * 1.0e7;
	}
	T energy = 0.0;
	T d_i = 1.0;
	T e_mass = 1.0;
	string flag = "membrane";
	ElectronTransmissionCoefficient
             (energy, transmission, conduction_band, d_i, e_mass,flag,_columnList);
	*/
	
	outputColumn(_columnList);
      }
    }
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
    for (int n=0; n<numMeshes; n++)
    {
      
        const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	//update charge field
	VectorT2Array& charge =
          dynamic_cast<VectorT2Array&>(_electricFields.charge[cells]);
        VectorT2Array& chargeN1 =
          dynamic_cast<VectorT2Array&>(_electricFields.chargeN1[cells]);

        if (_options.timeDiscretizationOrder > 1)
        {
            VectorT2Array& chargeN2 =
              dynamic_cast<VectorT2Array&>(_electricFields.chargeN2[cells]);
            chargeN2 = chargeN1;
        }
        chargeN1 = charge;
      	
	//update source term in electrostatics
	TArray& totalcharge = dynamic_cast<TArray&>(_electricFields.total_charge[cells]);
	for (int c=0; c<nCells; c++){
	  totalcharge[c] = - (charge[c][0] + charge[c][1]) * QE;
	  //totalcharge[c] = 0.0;
	}
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

	const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
	
	//if (vc.vcType == "dielectric"){
        
	  const StorageSite& cells = mesh.getCells();
	
	  MultiField::ArrayIndex cIndex(&_electricFields.charge,&cells);

	  ls.getX().addArray(cIndex,_electricFields.charge.getArrayPtr(cells));

	  const CRConnectivity& cellCells = mesh.getCellCells();
        
	  shared_ptr<Matrix> m(new CRMatrix<Tensor2x2,Tensor2x2,VectorT2>(cellCells));

	  ls.getMatrix().addMatrix(cIndex,cIndex,m);

	  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;

	      MultiField::ArrayIndex fIndex(&_electricFields.chargeFlux,&faces);
	      ls.getX().addArray(fIndex,_electricFields.chargeFlux.getArrayPtr(faces));

	      const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	      shared_ptr<Matrix> mft(new FluxJacobianMatrix<Tensor2x2,VectorT2>(faceCells));
	      ls.getMatrix().addMatrix(fIndex,cIndex,mft);

	      shared_ptr<Matrix> mff(new DiagonalMatrix<Tensor2x2,VectorT2>(faces.getCount()));
	      ls.getMatrix().addMatrix(fIndex,fIndex,mff);
	    }

	  foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;

	      MultiField::ArrayIndex fIndex(&_electricFields.chargeFlux,&faces);
	      ls.getX().addArray(fIndex,_electricFields.chargeFlux.getArrayPtr(faces));
	      
	      const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	      shared_ptr<Matrix> mft(new FluxJacobianMatrix<Tensor2x2,VectorT2>(faceCells));
	      ls.getMatrix().addMatrix(fIndex,cIndex,mft);

	      shared_ptr<Matrix> mff(new DiagonalMatrix<Tensor2x2,VectorT2>(faces.getCount()));
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
                                            _electricFields.potential_gradient));
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

            const ElectricBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _electricFields.potential,
                                  _electricFields.potential_flux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedPotential")
            {
                const T bT(bc["specifiedPotential"]);
                gbc.applyDirichletBC(bT);
            }
            else if (bc.bcType == "SpecifiedPotentialFlux")
            {
                const T specifiedFlux(bc["specifiedPotentialFlux"]);
                gbc.applyNeumannBC(specifiedFlux);
            }
	    else if ((bc.bcType == "Symmetry"))
            {
                T zeroFlux(NumTypeTraits<T>::getZero());
                gbc.applyNeumannBC(zeroFlux);
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
  
    if (_options.injection_enable){
      shared_ptr<Discretization>
	inj(new InjectionDiscretization<VectorT2, Tensor2x2, Tensor2x2>
	    (_meshes, _geomFields,
	     _electricFields.charge,
	     _electricFields.electric_field,
	     _electricFields.conduction_band,
	     _constants
	     ));
      discretizations.push_back(inj);
    }

  
     if (_options.tunneling_enable){
      shared_ptr<Discretization>
	tnd(new TunnelingDiscretization<VectorT2, Tensor2x2, Tensor2x2>
	    (_meshes, _geomFields,
	     _electricFields.charge,
	     _electricFields.electron_totaltraps,
	     _electricFields.conduction_band,
	     _constants
	     ));
      discretizations.push_back(tnd);      
    }


    if (_options.emission_enable){
      shared_ptr<Discretization>
	em(new EmissionDiscretization<VectorT2, Tensor2x2, Tensor2x2>
	   (_meshes, _geomFields,
	    _electricFields.charge,
	    _electricFields.electric_field,
	    _constants));
      discretizations.push_back(em);
    }

    if (_options.capture_enable){
      shared_ptr<Discretization>
	capt(new CaptureDiscretization<VectorT2, Tensor2x2, Tensor2x2>
	   (_meshes, _geomFields,
	    _electricFields.charge,
	    _electricFields.electron_totaltraps,
	    _electricFields.free_electron_capture_cross,
	    _constants));
      discretizations.push_back(capt);
    }
   
 
    if (_options.trapbandtunneling_enable){
      shared_ptr<Discretization>
	tbt(new TrapBandTunnelingDiscretization<VectorT2, Tensor2x2, Tensor2x2>
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
	dd(new ElecDiffusionDiscretization<VectorT2, Tensor2x2,Tensor2x2>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.diffusivity,
	    _electricFields.chargeGradient));
      discretizations.push_back(dd);
    }
#endif
    if (_options.drift_enable){
      shared_ptr<Discretization>
	cd(new DriftDiscretization<VectorT2, Tensor2x2, Tensor2x2>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.convectionFlux));
      discretizations.push_back(cd);
    }

    if (_options.transient_enable)
      {
      shared_ptr<Discretization>
	td(new TimeDerivativeDiscretization<VectorT2, Tensor2x2,Tensor2x2>
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
            

            GenericBCS<VectorT2,Tensor2x2,Tensor2x2> gbc(faces,mesh,
                                  _geomFields,
                                  _electricFields.charge,
                                  _electricFields.chargeFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());
	    //dielectric charging uses fixed zero dirichlet bc
	    VectorT2 zero;
	    zero[0] = zero[1] = T(0);
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
            GenericBCS<VectorT2,Tensor2x2,Tensor2x2> gbc(faces,mesh,
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
	    
	  MFRPtr eNormRatio((*eNorm)/(*_initialElectroStaticsNorm));

	  if (_options.printNormalizedResiduals)
	    cout << _niters << ": " << *eNormRatio << ";" <<  endl;
	  else
	    cout << _niters << ": " << *eNorm << ";" <<  endl;

	  if (*eNormRatio < _options.electrostaticsTolerance ) 
	    {
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
      
    if (flag1 && flag2 ) return true;   
    return false;
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
      const int nCells = cells.getCount();
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
      const int nCells = cells.getCount();
      const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array& > (_electricFields.electric_field[cells]);
      VectorT3Array& electron_velocity = dynamic_cast<VectorT3Array& > (_electricFields.electron_velocity[cells]);
      const T& electron_mobility = _constants["electron_mobility"];
      const T& electron_saturation_velocity =  _constants["electron_saturation_velocity"];

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
      
      const int nCells = cells.getCount();
      
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

    const T& membrane_workfunction = _constants["membrane_workfunction"];
    const T& substrate_workfunction = _constants["substrate_workfunction"];
    const T& dielectric_thickness = _constants["dielectric_thickness"];
    const T& optical_dielectric_constant = _constants["optical_dielectric_constant"];
    const T& dielectric_ionization = _constants["dielectric_ionization"];
    const T& electron_trapdepth = _constants["electron_trapdepth"];
    const T& temperature = _constants["OP_temperature"];
    const T& electron_effmass = _constants["electron_effmass"];
    const T& poole_frenkel_emission_frequency = _constants["poole_frenkel_emission_frequency"];

    //??? dielectric thickness ///
    T effefield = (membrane_workfunction - substrate_workfunction) / dielectric_thickness;

    T alpha = sqrt(QE / (PI * E0_SI * optical_dielectric_constant));

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      const ElectricVC<T>& vc = *_vcMap[mesh.getID()];
      if (vc.vcType == "dielectric"){

	const StorageSite& cells = mesh.getCells(); 
	const int nCells = cells.getCount();

	const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);
	const TArray& electron_totaltraps = dynamic_cast<TArray& > (_electricFields.electron_totaltraps[cells]);
	TArray& free_electron_capture_cross = dynamic_cast<TArray&> (_electricFields.free_electron_capture_cross[cells]);
	//  charge is a VectorT2 array
	//  charge[0] contains the electron density in traps
	//  charge[1] contains the electron density in conduction band 

	VectorT2Array& charge = dynamic_cast<VectorT2Array& > (_electricFields.charge[cells]);
 	VectorT2Array& chargeN1 = dynamic_cast<VectorT2Array& > (_electricFields.chargeN1[cells]);
        VectorT2Array* chargeN2  = 0;
	if (_options.timeDiscretizationOrder > 1){
	  chargeN2 = &(dynamic_cast<VectorT2Array& > (_electricFields.chargeN2[cells]));
	}

	for(int c=0; c<nCells; c++)
	{ 	  
	  // need to figure out what direction to use ???//
	  T fermilevel = -substrate_workfunction+effefield*cellCentroid[c][2];  
	  T energy(0.);

	  //fermi function is extremely sentitive to energy and fermilevel ???//

	  energy = -dielectric_ionization - electron_trapdepth;
	  
	  charge[c][0] = electron_totaltraps[c] * FermiFunction(energy, fermilevel, temperature);
	  chargeN1[c][0] = charge[c][0];

	  //here it should be N2 but got a compile error
	  if (_options.timeDiscretizationOrder > 1)
	    (*chargeN2)[c][0] = charge[c][0];

	  energy = -dielectric_ionization;
	  
	  charge[c][1] = electron_totaltraps[c] * FermiFunction(energy, fermilevel, temperature);
	  chargeN1[c][1] = charge[c][1];
	 
	  
	  if (_options.timeDiscretizationOrder > 1)
	    (*chargeN2)[c][1] = charge[c][1];

	  T expt = (electron_trapdepth - alpha * sqrt(fabs(effefield))) * QE / (K_SI*temperature);
	  T beta = exp(-expt);
	  T velocity = sqrt(8 * K_SI * temperature / (PI * ME * electron_effmass));

	  free_electron_capture_cross[c] = charge[c][0] * poole_frenkel_emission_frequency * beta;
	  free_electron_capture_cross[c] /= (velocity*(electron_totaltraps[c]-charge[c][0])*charge[c][1]);
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
    
        ibP->zero();

	mIC.multiplyAndAdd(*ibP,cP);
	mIP.multiplyAndAdd(*ibP,sP);
#if 0
	const Array<int>& ibFaceList = mesh.getIBFaceList();
	const StorageSite& faces = mesh.getFaces();
	
	for(int f=0; f<ibFaces.getCount();f++){
	  int fID = ibFaceList[f];
	  (*ibP)[fID] = 10.0;
	}
#endif
        _electricFields.potential.addArray(ibFaces,ibP);
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

    const VectorT3Array& solidFaceCoordinate =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[solidFaces]);

       
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
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
	for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)	    {
	  const int c = sFCCol[nc];
	  const T coeff = iCoeffs[nc];
	  const T efmag = mag(electric_field[c]);
	  forceSign = dot(electric_field[c], solidFaceArea[f]);
	  if (fabs(forceSign) > 0.0) 
	    forceSign /= fabs(forceSign);
	  else{
	    forceSign = 0.0;
	    cout << f << " " << solidFaceCoordinate[f][0] << " "  <<solidFaceCoordinate[f][1] << " " <<efmag << endl;
	  }
	  forceMag += 0.5 * coeff * dielectric_constant[c] *  efmag * efmag * forceSign; 
	}
	const VectorT3& Af = solidFaceArea[f];
	force[f] = Af * forceMag;
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




    


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  ElectricFields& _electricFields;

  ElectricBCMap _bcMap;
  ElectricVCMap _vcMap;

  ElectricModelOptions<T> _options; 
  ElectricModelConstants<T> _constants;
  GradientModel<T> _potentialGradientModel;
  GradientModel<VectorT2> _chargeGradientModel;
  
  MFRPtr _initialElectroStaticsNorm;
  MFRPtr _initialChargeTransportNorm;
  int _niters;

  StorageSite _columns;
  ColumnList _columnList;

  shared_ptr<CRConnectivity> _cellColumns;
  
  shared_ptr<CRConnectivity> _columnCells;
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
  _impl->advance(niter);
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
void
ElectricModel<T>::computeSolidSurfaceForcePerUnitArea(const StorageSite& particles)
{
  return _impl->computeSolidSurfaceForce(particles,true);
}
