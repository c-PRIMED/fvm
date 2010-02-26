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
#include "ConvectionDiscretization.h"
#include "TimeDerivativeDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "SourceDiscretization.h"
#include "Octree.h"

template<class T>
class ElectricModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  
  typedef Gradient<T> CGradType;
  typedef Array<Gradient<T> > CGradArray;

  Impl(const GeomFields& geomFields,
       ElectricFields& electricFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _electricFields(electricFields),
    _potentialGradientModel(_meshes,_electricFields.potential,
                              _electricFields.electric_field,_geomFields),
    _chargeGradientModel (_meshes, _electricFields.charge,
			  _electricFields.chargeGradient,_geomFields),
    _initialElectroStaticsNorm(),
    _initialChargeTransportNorm(),
    _initialTunnelingTransportNorm(),
    _niters(0)
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
	
	if (_options.electrostatics_enable){
	
	  //initial potential setup
	  shared_ptr<TArray> pCell(new TArray(cells.getCount()));
	  *pCell = _options["initialPotential"];
	  _electricFields.potential.addArray(cells,pCell);

	  //dielectric_constant setup
	  shared_ptr<TArray> permCell(new TArray(cells.getCount()));
	  *permCell = vc["dielectric_constant"];
          _electricFields.dielectric_constant.addArray(cells,permCell);
	  
	  //total_charge (source in Poisson equation) setup
	  shared_ptr<TArray> sCell(new TArray(cells.getCount()));
	  *sCell = _options["initialTotalCharge"];
	  _electricFields.total_charge.addArray(cells,sCell);
	
	  //electric_field (potential gradient vector) setup
	  shared_ptr<CGradArray> gradp(new CGradArray(cells.getCount()));
          gradp->zero();	
          _electricFields.electric_field.addArray(cells,gradp);

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

	if (_options.chargetransport){

	  //initial charge density
	  shared_ptr<TArray> cCell(new TArray(cells.getCount()));
	  *cCell = _options["initialCharge"];
	  _electricFields.charge.addArray(cells,cCell);
	
	  if (_options.transient)
	    {
	      _electricFields.chargeN1.addArray(cells,
						  dynamic_pointer_cast<ArrayBase>(cCell->newCopy()));
	      if (_options.timeDiscretizationOrder > 1)
		_electricFields.chargeN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(cCell->newCopy()));
	    }

	  //diffusivity 
	  shared_ptr<TArray> diffCell(new TArray(cells.getCount()));
	  *diffCell = 1.0;
	  _electricFields.diffusivity.addArray(cells,diffCell);
	
	//create a zero field
	  shared_ptr<TArray> zeroCell(new TArray(cells.getCount()));
	  *zeroCell = 0.0;
	  _electricFields.zero.addArray(cells,zeroCell);

	//create a one field
	  shared_ptr<TArray> oneCell(new TArray(cells.getCount()));
	  *oneCell = 1.0;
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
	      
	      shared_ptr<TArray> cFlux(new TArray(faces.getCount()));
	      cFlux->zero();
	      _electricFields.chargeFlux.addArray(faces,cFlux);
          
	    }
	  foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<TArray> cFlux(new TArray(faces.getCount()));
	      cFlux->zero();
	      _electricFields.chargeFlux.addArray(faces,cFlux);
          
	    }
	
	  //inital convection flux at faces
	  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<TArray> convFlux(new TArray(faces.getCount()));
	      convFlux->zero();
	      _electricFields.convectionFlux.addArray(faces,convFlux);
          
	    }
	  foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	    {
	      const FaceGroup& fg = *fgPtr;
	      const StorageSite& faces = fg.site;
	      
	      shared_ptr<TArray> convFlux(new TArray(faces.getCount()));
	      convFlux->zero();
	      _electricFields.convectionFlux.addArray(faces,convFlux);
	      
	    }

	}
	
	if (_options.tunneling){
	  shared_ptr<TArray> tcCell(new TArray(cells.getCount()));
	  *tcCell = _options["initialTunnelingCharge"];
	  _electricFields.tunnelingCharge.addArray(cells,tcCell);
	
	  
	}
    }

    _niters  = 0;
    _initialElectroStaticsNorm = MFRPtr();
    _initialChargeTransportNorm = MFRPtr();
    _initialTunnelingTransportNorm = MFRPtr();
  }
  
  ElectricBCMap& getBCMap() {return _bcMap;}

  ElectricBC<T>& getBC(const int id) {return *_bcMap[id];}

  ElectricVCMap& getVCMap() {return _vcMap;}

  ElectricVC<T>& getVC(const int id) {return *_vcMap[id];}
  
  ElectricModelOptions<T>& getOptions() {return _options;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      
        const Mesh& mesh = *_meshes[n];
	//update charge field
        const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	/*
	TArray& charge =
          dynamic_cast<TArray&>(_electricFields.charge[cells]);
        TArray& chargeN1 =
          dynamic_cast<TArray&>(_electricFields.chargeN1[cells]);

        if (_options.timeDiscretizationOrder > 1)
        {
            TArray& chargeN2 =
              dynamic_cast<TArray&>(_electricFields.chargeN2[cells]);
            chargeN2 = chargeN1;
        }
        chargeN1 = charge;
      	
	//update tunneling field
	TArray& tunneling = dynamic_cast<TArray&>(_electricFields.tunnelingCharge[cells]);
      
	//update source term in electrostatics
	TArray& source = dynamic_cast<TArray&>(_electricFields.source[cells]);
	for (int c=0; c<nCells; c++){
	  source[c]=tunneling[c];
	}
	*/
        
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

  MFRPtr solveTunnelingTransport()
  {

   //calculate tunneling charge
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
      const Mesh& mesh = *_meshes[n];
      
      const StorageSite& cells = mesh.getCells();

      const int nCells = cells.getCount();

      TArray& tunnelingCharge = dynamic_cast<TArray&>(_electricFields.tunnelingCharge[cells]);
      const TArray& potential = dynamic_cast<const TArray&>(_electricFields.potential[cells]);
      const VectorT3Array& elecField = dynamic_cast<const VectorT3Array&>(_electricFields.electric_field[cells]);

      for(int c=0; c<nCells; c++){
  	tunnelingCharge[c]=1.0;
      }
    }
    
    //return the residual of tunneling charge field
    MFRPtr tNorm;
    
    return tNorm;

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

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex cIndex(&_electricFields.charge,&cells);

        ls.getX().addArray(cIndex,_electricFields.charge.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(cIndex,cIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_electricFields.chargeFlux,&faces);
            ls.getX().addArray(fIndex,_electricFields.chargeFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,cIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_electricFields.chargeFlux,&faces);
            ls.getX().addArray(fIndex,_electricFields.chargeFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,cIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,cIndex,mff);
        }

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
                                            _electricFields.electric_field));
    discretizations.push_back(dd);    

    shared_ptr<Discretization>
      sd(new SourceDiscretization<T>
	 (_meshes, 
	  _geomFields, 
	  _electricFields.potential,
	  _electricFields.total_charge
	));
    discretizations.push_back(sd);

    if(_options.ibm){
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
    if (_options.diffusion){
      shared_ptr<Discretization>
	dd(new DiffusionDiscretization<T,T,T>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.diffusivity,
	    _electricFields.chargeGradient));
      discretizations.push_back(dd);
    }
    if (_options.drift){
      shared_ptr<Discretization>
	cd(new ConvectionDiscretization<T,T,T>
	   (_meshes,_geomFields,
	    _electricFields.charge,
	    _electricFields.convectionFlux,
	    _electricFields.zero,
	    _electricFields.chargeGradient));
      discretizations.push_back(cd);
    }
    if (_options.transient)
    {
      shared_ptr<Discretization>
	td(new TimeDerivativeDiscretization<T,T,T>
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
                                  _electricFields.charge,
                                  _electricFields.chargeFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            if (bc.bcType == "SpecifiedCharge")
            {
                const T bT(bc["specifiedCharge"]);
                gbc.applyDirichletBC(bT);
            }

            else if (bc.bcType == "SpecifiedChargeFlux")
            {
	       const T specifiedFlux(bc["specifiedChargeFlux"]);
	       gbc.applyNeumannBC(specifiedFlux);
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
                                  _electricFields.charge,
                                  _electricFields.chargeFlux,
                                  ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }
    }
  }



 
  bool advance(const int niter)
  {
    for(int n=0; n<niter; n++)
    {
      bool flag1 = false;
      bool flag2 = false;
      bool flag3 = false; 
      
      if(_options.electrostatics_enable){

        MFRPtr eNorm = solveElectroStatics();
	if (_niters < 5)
        {
	  _initialElectroStaticsNorm->setMax(*eNorm);
	}
	MFRPtr eNormRatio((*eNorm)/(*_initialElectroStaticsNorm));

	if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *eNormRatio << ";" <<  endl;
        else
	  cout << _niters << ": " << *eNorm << ";" <<  endl;

	if (*eNormRatio < _options.electrostaticsTolerance) 
	  flag1 = true;
      }

       
      if(_options.chargetransport){
     
	MFRPtr cNorm = solveChargeTransport();
        if (_niters < 5)
        {
	  _initialChargeTransportNorm->setMax(*cNorm);
        }
               
	MFRPtr cNormRatio((*cNorm)/(*_initialChargeTransportNorm));
        
        if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *cNormRatio <<  endl;
        else
	  cout << _niters << ": " <<  *cNorm <<  endl;

	if (*cNormRatio < _options.chargetransportTolerance)  
	  flag2 = true;
      }

      if(_options.tunneling){
	/*

	MFRPtr tNorm = solveTunnelingTransport();
	if (_niters < 5){
	  _initialTunnelingTransportNorm->setMax(*tNorm);
	}
	MFRPtr tNormRatio((*tNorm)/(*_initialTunnelingTransportNorm));
        
        if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *tNormRatio <<  endl;
        else
	  cout << _niters << ": " <<  *tNorm <<  endl;

	if (*tNormRatio < _options.tunnelingtransportTolerance)
	  flag3=true;
	*/
      }
	
      _niters++;
      
      if (flag1 && flag2 && flag3) return true;    
     
    }
    return false;
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


 const int findClosestPoint(const VectorT3 point, Octree& O)
 {
   const int nearestCell = O.getNode(point);
   return nearestCell;
  
 }

 vector<int> findCloestPoints(VectorT3 point, Octree& O, const double radius)
 {
   vector<int> nearestCellList;
   O.getNodes(point, radius, nearestCellList);
 }

 shared_ptr<VectorT3Array> createPathAndDiscretize (
	 const VectorT3 startPoint, const VectorT3 endPoint, const int N)
 {
   shared_ptr<VectorT3Array> pathPoints(new VectorT3Array(N+1));
   const VectorT3 path = endPoint - startPoint;
   const VectorT3 delta(path / T(N));
   for (int n=0; n<=N; n++){
     VectorT3 tmp = startPoint + (delta * T(n));
     (*pathPoints)[n] = tmp;
   }
   return pathPoints;
 }


 shared_ptr<VectorT3Array> pathPointsInterpolation(
			      shared_ptr<VectorT3Array> pathPoints, 
			      Octree& O,
			      const double radius)
 {
   const int  printOption = 1;
   const int nP = (*pathPoints).getLength();
   shared_ptr<VectorT3Array> pathPointValues (new VectorT3Array(nP));
   for(int n=0; n<nP; n++){
     (*pathPointValues)[n]=0.0;
   }

   const Mesh& mesh = *_meshes[0];
   const StorageSite& cells = mesh.getCells();
   const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);

   /// print out path Points 
   if (printOption == 1) {
     for (int n=0; n<nP; n++){
       cout<<n<<"   "<<(*pathPoints)[n]<<endl;
     }
     cout<<"end of path points list"<<endl;
   }
      
   for (int n=0; n<nP; n++){ 

     //step 1. for each path point, find out its cell neighbors within a radius
     const VectorT3 point = (*pathPoints)[n];
     vector<int> neighborList;
     O.getNodes(point, radius, neighborList);

     //step 2. use distance weighted method, calcualte the interpolation coefficients to all the neighbors
     const int nNBCells = neighborList.size();
     TArray pointToCellCoeff(nNBCells);
     T wtSum(0);
     int onTarget = 0;
     int target = 0;
     for (int i=0; i<nNBCells; i++){
       const int cid = neighborList[i];
       cout<<"point "<<n<<"   "<<point<<endl;
       cout<<"neighbors "<<i<<"   "<<cellCentroid[cid]<<endl;

       VectorT3 dr(point - cellCentroid[cid]);
       if (dot(dr,dr)<=1.0e-15){
	 onTarget = 1;
	 target = i;
       }
       else{
	 T wt = 1.0/dot(dr,dr);
	 pointToCellCoeff[i]=wt;
	 wtSum += wt;
       }
     }
     if (onTarget == 1){
       for (int i=0; i<nNBCells; i++){
	 pointToCellCoeff[i] = 0;
       }
       pointToCellCoeff[target] = 1;
     } 
     else{
       for (int i=0; i<nNBCells; i++){
	 pointToCellCoeff[i] /= wtSum;
       }
     }

     //step 3. interpolation the value using the coefficients
     for (int i=0; i<nNBCells; i++){
       const int cid = neighborList[i];
       VectorT3 cellValue;
       cellValue[0]=1.0;
       cellValue[1]=1.0;
       cellValue[2]=1.0;       
       (*pathPointValues)[n] = (*pathPointValues)[n] + cellValue * pointToCellCoeff[i];
       cout<<"coeff   "<<pointToCellCoeff[i]<<endl;
     } 
   }

    /// print out path point values 
   if (printOption){
     for (int n=0; n<nP; n++){
       cout<<n<<"   "<<(*pathPointValues)[n]<<endl;
     }
     cout<<"end of path points list"<<endl;
   }
   return pathPointValues;
 }



   
   


    
     


private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  ElectricFields& _electricFields;

  ElectricBCMap _bcMap;
  ElectricVCMap _vcMap;

  ElectricModelOptions<T> _options; 
  GradientModel<T> _potentialGradientModel;
  GradientModel<T> _chargeGradientModel;
  
  MFRPtr _initialElectroStaticsNorm;
  MFRPtr _initialChargeTransportNorm;
  MFRPtr _initialTunnelingTransportNorm;
  int _niters;
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
#if 0
template<class T>
void
ElectricModel<T>::computeIBFacePotential(const StorageSite& particles)
{
  _impl->computeIBFacePotential(particles);
}
#endif 

/* only the members who need to be exposed in python are listed here */
/* and declared in ElectricModel class */

/*

template<class T>
const int
ElectricModel<T>::findClosestPoint(const VectorT3 point, Octree& O)
{
  return _impl->findClosestPoint(point, O);
}
 
template<class T>
boost::shared_ptr<ArrayBase>
ElectricModel<T>::createPathAndDiscretize 
    (const VectorT3 startPoint, const VectorT3 endPoint, const int N)
{
  return _impl->createPathAndDiscretize(startPoint, endPoint, N);
}

template<class T>
boost::shared_ptr<ArrayBase>
ElectricModel<T>::pathPointsInterpolation(boost::shared_ptr<ArrayBase> pathPoints, Octree& O, const double radius)
{
  return _impl->pathPointsInterpolation( dynamic_pointer_cast<Array<Vector<T,3> > >(pathPoints), O, radius);
}
*/
