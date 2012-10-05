// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include "Mesh.h"

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>

#include "NumType.h"
#include "Array.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"
#include "GenericBCS.h"
#include "Vector.h"
#include "VectorTranspose.h"
#include "DiffusionDiscretization.h"
#include "ConvectionDiscretization.h"
#include "TimeDerivativeDiscretization.h"
#include "Underrelaxer.h"
#include "MomentumPressureGradientDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "GenericIBDiscretization.h"
#include "StressTensor.h"

template<class T>
class FlowModel<T>::Impl
{
public:

  typedef Array<int> IntArray;
  
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef VectorTranspose<T,3> VectorT3T;

  typedef Array<VectorT3> VectorT3Array;
  typedef DiagonalTensor<T,3> DiagTensorT3;

  typedef CRMatrix<DiagTensorT3,T,VectorT3> VVMatrix;
  typedef typename VVMatrix::DiagArray VVDiagArray;
  
  typedef CRMatrix<T,T,T> PPMatrix;
  typedef typename PPMatrix::DiagArray PPDiagArray;
  typedef typename PPMatrix::PairWiseAssembler PPAssembler;

#ifdef PV_COUPLED
  typedef CRMatrixRect<VectorT3T,VectorT3,T> PVMatrix;
  typedef typename PVMatrix::DiagArray PVDiagArray;
  typedef typename PVMatrix::PairWiseAssembler PVAssembler;

  typedef CRMatrixRect<VectorT3,T,VectorT3> VPMatrix;
  typedef typename VPMatrix::DiagArray VPDiagArray;
  typedef typename VPMatrix::PairWiseAssembler VPAssembler;

#endif
  
  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;

  typedef Gradient<T> PGradType;
  typedef Array<PGradType> PGradArray;

  typedef FluxJacobianMatrix<T,T> FMatrix;
  
  typedef StressTensor<T> StressTensorT6;
  typedef Array<StressTensorT6> StressTensorArray;
  
  Impl(const GeomFields& geomFields,
       FlowFields& thermalFields,
       const MeshList& meshes):
    _meshes(meshes),
    _geomFields(geomFields),
    _flowFields(thermalFields),
    _velocityGradientModel(_meshes,_flowFields.velocity,
                           _flowFields.velocityGradient,_geomFields),
    _pressureGradientModel(_meshes,_flowFields.pressure,
                           _flowFields.pressureGradient,_geomFields),
    _initialMomentumNorm(),
    _initialContinuityNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        FlowVC<T> *vc(new FlowVC<T>());
        vc->vcType = "flow";
        _vcMap[mesh.getID()] = vc;
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            if (_bcMap.find(fg.id) == _bcMap.end())
            {
                FlowBC<T> *bc(new FlowBC<T>());
                
                _bcMap[fg.id] = bc;
                if ((fg.groupType == "wall"))
                {
                    bc->bcType = "NoSlipWall";
                }
                else if ((fg.groupType == "velocity-inlet"))
                {
                    bc->bcType = "VelocityBoundary";
                }
                else if ((fg.groupType == "pressure-inlet") ||
                         (fg.groupType == "pressure-outlet"))
                {
                    bc->bcType = "PressureBoundary";
                }
                else if ((fg.groupType == "symmetry"))
                {
                    bc->bcType = "Symmetry";
                }
		else
                  throw CException("FlowModel: unknown face group type "
                                   + fg.groupType);
            }
        }
	/*
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
	    if ((fg.groupType == "ESinterface"))
	      { bc->bcType = "ESInterfaceBC";}
	  }
	*/
	

    }
  }

  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const FlowVC<T>& vc = *_vcMap[mesh.getID()];
            
        const StorageSite& cells = mesh.getCells();
        const StorageSite& faces = mesh.getFaces();

        shared_ptr<VectorT3Array> vCell(new VectorT3Array(cells.getCountLevel1()));

        VectorT3 initialVelocity;
        initialVelocity[0] = _options["initialXVelocity"];
        initialVelocity[1] = _options["initialYVelocity"];
        initialVelocity[2] = _options["initialZVelocity"];
        *vCell = initialVelocity;

        _flowFields.velocity.addArray(cells,vCell);

        shared_ptr<VectorT3Array> uparallelCell(new VectorT3Array(cells.getCount()));
        uparallelCell ->zero();
        _flowFields.uparallel.addArray(cells,uparallelCell);

      shared_ptr<VectorT3Array> tauCell(new VectorT3Array(faces.getCount()));
        tauCell ->zero();
        _flowFields.tau.addArray(faces,tauCell);

     shared_ptr<VectorT3Array> TauWallCell(new VectorT3Array(faces.getCount()));
        TauWallCell ->zero();
        _flowFields.tauwall.addArray(faces,TauWallCell);




  
       

        if (_options.transient)
        {
            _flowFields.velocityN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(vCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
              _flowFields.velocityN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(vCell->newCopy()));

        }
        
        shared_ptr<TArray> pCell(new TArray(cells.getCountLevel1()));
        shared_ptr<TArray> pFace(new TArray(faces.getCountLevel1()));
        *pCell = _options["initialPressure"];
        *pFace = _options["initialPressure"];
        _flowFields.pressure.addArray(cells,pCell);
        _flowFields.pressure.addArray(faces,pFace);


        shared_ptr<TArray> rhoCell(new TArray(cells.getCountLevel1()));
        *rhoCell = vc["density"];
        _flowFields.density.addArray(cells,rhoCell);

        shared_ptr<TArray> muCell(new TArray(cells.getCountLevel1()));
        *muCell = vc["viscosity"];
        _flowFields.viscosity.addArray(cells,muCell);

        shared_ptr<PGradArray> gradp(new PGradArray(cells.getCountLevel1()));
        gradp->zero();
        _flowFields.pressureGradient.addArray(cells,gradp);

        shared_ptr<TArray> ci(new TArray(cells.getCountLevel1()));
        ci->zero();
        _flowFields.continuityResidual.addArray(cells,ci);

        // compute default value of mass flux
        
        const VectorT3Array& faceArea =
          dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
        const VectorT3Array& V = *vCell;
        const TArray& density = *rhoCell;
        
        const CRConnectivity& faceCells = mesh.getAllFaceCells();

        const int nFaces = faces.getCount();

        shared_ptr<TArray> mfPtr(new TArray(faces.getCount()));
        TArray& mf = *mfPtr;

        for(int f=0; f<nFaces; f++)
        {
            const int c0 = faceCells(f,0);
            const int c1 = faceCells(f,1);
            
            mf[f] = 0.5*(density[c0]*dot(V[c0],faceArea[f]) +
                         density[c1]*dot(V[c1],faceArea[f]));
        }
        _flowFields.massFlux.addArray(faces,mfPtr);

        // store momentum flux at interfaces
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            shared_ptr<VectorT3Array> momFlux(new VectorT3Array(faces.getCount()));
            momFlux->zero();
            _flowFields.momentumFlux.addArray(faces,momFlux);
	    
	    if(fg.groupType == "ESinterface"){
	      cout << "interface init" <<endl;
	      const StorageSite& Intfaces = fg.site; 
	      
	      shared_ptr<VectorT3Array> InterfaceVelFace(new VectorT3Array(Intfaces.getCount()));
	      InterfaceVelFace ->zero();
	      _flowFields.InterfaceVelocity.addArray(Intfaces,InterfaceVelFace);
	      
	      shared_ptr<StressTensorArray> InterfaceStressFace(new StressTensorArray(Intfaces.getCount()));
	      InterfaceStressFace ->zero();
	      _flowFields.InterfaceStress.addArray(Intfaces,InterfaceStressFace);
	      
	      shared_ptr<TArray> InterfacePressFace(new TArray(Intfaces.getCount()));
	      *InterfacePressFace = _options["initialPressure"];
	      _flowFields.InterfacePressure.addArray(Intfaces,InterfacePressFace);
	      
	      shared_ptr<TArray> InterfaceDensityFace(new TArray(Intfaces.getCount()));
	      *InterfaceDensityFace =vc["density"];
	      _flowFields.InterfaceDensity.addArray(Intfaces,InterfaceDensityFace);
	      
	      
	    }
	    
	    
        }
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            const FlowBC<T>& bc = *_bcMap[fg.id];

            FloatValEvaluator<VectorT3>
              bVelocity(bc.getVal("specifiedXVelocity"),
                        bc.getVal("specifiedYVelocity"),
                        bc.getVal("specifiedZVelocity"),
                        faces);

            const int nFaces = faces.getCount();
            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            if ((bc.bcType == "NoSlipWall") ||
                (bc.bcType == "SlipJump") ||
                (bc.bcType == "VelocityBoundary"))
            {
                // arrays for this face group
                TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
                
                const VectorT3Array& faceArea =
                  dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
                
                for(int f=0; f<nFaces; f++)
                {
                    const int c0 = faceCells(f,0);
                    massFlux[f] = density[c0]*dot(bVelocity[f],faceArea[f]);
                 }
            }
            else if (bc.bcType == "SpecifiedPressure")
            {
                T bp = bc["specifiedPressure"];
                TArray& facePressure = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
                TArray& cellPressure = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
                for(int f=0; f<nFaces; f++)
                {
                    const int c1 = faceCells(f,1);
                    facePressure[f] = cellPressure[c1] = bp;
                }
            }
            else if ((bc.bcType == "Symmetry"))
            {
                TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
                massFlux.zero();
            }
            shared_ptr<VectorT3Array> momFlux(new VectorT3Array(faces.getCount()));
            momFlux->zero();
            _flowFields.momentumFlux.addArray(faces,momFlux);
           
	    
	    
        }
    }

    computeContinuityResidual();
    _niters  =0;
    _initialMomentumNorm = MFRPtr();
    _initialContinuityNorm = MFRPtr();
  }
  
  FlowBCMap& getBCMap() {return _bcMap;}
  FlowVCMap& getVCMap() {return _vcMap;}
  FlowModelOptions<T>& getOptions() {return _options;}

  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        VectorT3Array& v =
          dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
        VectorT3Array& vN1 =
          dynamic_cast<VectorT3Array&>(_flowFields.velocityN1[cells]);

        if (_options.timeDiscretizationOrder > 1)
        {
            VectorT3Array& vN2 =
              dynamic_cast<VectorT3Array&>(_flowFields.velocityN2[cells]);
            vN2 = vN1;
        }
        vN1 = v;
    }
  }
 const Field& getViscosityField() const
  {
    if (_options.turbulent)
      return _flowFields.totalviscosity;
    else
      return _flowFields.viscosity;
  }

  void computeIBFaceVelocity(const StorageSite& particles)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    typedef CRMatrixTranspose<T,VectorT3,VectorT3> IMatrixV3;

    const VectorT3Array& pV =
      dynamic_cast<const VectorT3Array&>(_flowFields.velocity[particles]);

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

           IMatrixV3 mICV(mIC);
	   

           GeomFields::SSPair key2(&ibFaces,&particles);
           const IMatrix& mIP =
           dynamic_cast<const IMatrix&>
            (*_geomFields._interpolationMatrices[key2]);

           IMatrixV3 mIPV(mIP);
	   

           shared_ptr<VectorT3Array> ibV(new VectorT3Array(ibFaces.getCount()));
        
           const VectorT3Array& cV =
            dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
    

           ibV->zero();

           mICV.multiplyAndAdd(*ibV,cV);
   	   mIPV.multiplyAndAdd(*ibV,pV);

#if 0
        ofstream debugFile;
	stringstream ss(stringstream::in | stringstream::out);
        ss <<  MPI::COMM_WORLD.Get_rank();
	string  fname1 = "IBVelocity_proc" +  ss.str() + ".dat";
	debugFile.open(fname1.c_str());
	
	//debug use
	const Array<int>& ibFaceList = mesh.getIBFaceList();
	const StorageSite& faces = mesh.getFaces();
	const VectorT3Array& faceCentroid = 
          dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[faces]);
	const double angV = 1.0;
	VectorT3 center;
	center[0]=0.;
	center[1]=0.;
	center[2]=0.;	

	for(int f=0; f<ibFaces.getCount();f++){
	  int fID = ibFaceList[f];
	  debugFile << "f=" <<   f << setw(10) <<  "   fID = " <<  fID << "  faceCentroid = " << faceCentroid[fID] << " ibV = " << (*ibV)[f] << endl;
	}
	  
	 debugFile.close();
#endif


          _flowFields.velocity.addArray(ibFaces,ibV);

	}
    }

  }

  map<string,shared_ptr<ArrayBase> >&
  getPersistenceData()
  {
    _persistenceData.clear();
    
    Array<int>* niterArray = new Array<int>(1);
    (*niterArray)[0] = _niters;
    _persistenceData["niters"]=shared_ptr<ArrayBase>(niterArray);
    
    if (_initialMomentumNorm)
    {
        _persistenceData["initialMomentumNorm"] =
          _initialMomentumNorm->getArrayPtr(_flowFields.velocity);
    }
    else
    {
        Array<Vector<T,3> >* xArray = new Array<Vector<T,3> >(1);
        xArray->zero();
        _persistenceData["initialMomentumNorm"]=shared_ptr<ArrayBase>(xArray);
        
    }
    
    if (_initialContinuityNorm)
    {
        _persistenceData["initialContinuityNorm"] =
          _initialContinuityNorm->getArrayPtr(_flowFields.pressure);
    }
    else
    {
        Array<T>* xArray = new Array<T>(1);
        xArray->zero();
        _persistenceData["initialContinuityNorm"]=shared_ptr<ArrayBase>(xArray);
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

    if (_persistenceData.find("initialMomentumNorm") != _persistenceData.end())
    {
        shared_ptr<ArrayBase>  r = _persistenceData["initialMomentumNorm"];
        _initialMomentumNorm = MFRPtr(new MultiFieldReduction());
        _initialMomentumNorm->addArray(_flowFields.velocity,r);
    }

    if (_persistenceData.find("initialContinuityNorm") != _persistenceData.end())
    {
        shared_ptr<ArrayBase>  r = _persistenceData["initialContinuityNorm"];
        _initialContinuityNorm = MFRPtr(new MultiFieldReduction());
        _initialContinuityNorm->addArray(_flowFields.pressure,r);
    }
    
    computeContinuityResidual();
    
  }
  

  
  void initMomentumLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);

        ls.getX().addArray(vIndex,_flowFields.velocity.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<DiagTensorT3,T,VectorT3>(cellCells));

        ls.getMatrix().addMatrix(vIndex,vIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_flowFields.momentumFlux,&faces);
            ls.getX().addArray(fIndex,_flowFields.momentumFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<DiagTensorT3,VectorT3>(faceCells));
            ls.getMatrix().addMatrix(fIndex,vIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<DiagTensorT3,VectorT3>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
        
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_flowFields.momentumFlux,&faces);
            ls.getX().addArray(fIndex,_flowFields.momentumFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<DiagTensorT3,VectorT3>(faceCells));
            ls.getMatrix().addMatrix(fIndex,vIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<DiagTensorT3,VectorT3>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
    }
  }

  void linearizeMomentum(LinearSystem& ls)
  {
    _velocityGradientModel.compute();
    
    DiscrList discretizations;
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<VectorT3,DiagTensorT3,T>
         (_meshes,_geomFields,
          _flowFields.velocity,
           getViscosityField(),
          _flowFields.velocityGradient));

    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<VectorT3,DiagTensorT3,T>
         (_meshes,_geomFields,
          _flowFields.velocity,
          _flowFields.massFlux,
          _flowFields.continuityResidual,
          _flowFields.velocityGradient));

    shared_ptr<Discretization>
      pd(new MomentumPressureGradientDiscretization<T>
         (_meshes,_geomFields,
          _flowFields,
          _pressureGradientModel));

    discretizations.push_back(dd);
    discretizations.push_back(cd);
    discretizations.push_back(pd);

    if (_options.transient)
    {
        shared_ptr<Discretization>
          td(new TimeDerivativeDiscretization<VectorT3,DiagTensorT3,T>
             (_meshes,_geomFields,
              _flowFields.velocity,
              _flowFields.velocityN1,
              _flowFields.velocityN2,
              _flowFields.density,
              _options["timeStep"]));
        
        discretizations.push_back(td);
    }


    shared_ptr<Discretization>
      ibm(new GenericIBDiscretization<VectorT3,DiagTensorT3,T>
             (_meshes,_geomFields,_flowFields.velocity));
      
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
            const int nFaces = faces.getCount();

            const FlowBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<VectorT3,DiagTensorT3,T> gbc(faces,mesh,
                                                    _geomFields,
                                                    _flowFields.velocity,
                                                    _flowFields.momentumFlux,
                                                    ls.getMatrix(),
                                                    ls.getX(),
                                                    ls.getB());

            const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);
            //            const CRConnectivity& faceCells = mesh.getFaceCells(faces);
    
            FloatValEvaluator<VectorT3>
              bVelocity(bc.getVal("specifiedXVelocity"),
                        bc.getVal("specifiedYVelocity"),
                        bc.getVal("specifiedZVelocity"),
                        faces);

            if (bc.bcType == "NoSlipWall")
            {
            	gbc.applyDirichletBC(bVelocity);
            }
            else if (bc.bcType == "SlipJump")
            {
                slipJumpMomentumBC(faces,mesh,
                                   gbc,
                                   bc["accomodationCoefficient"],
                                   bVelocity);
            }
            else if ((bc.bcType == "VelocityBoundary") ||
                     (bc.bcType == "PressureBoundary"))
            {
                for(int f=0; f<nFaces; f++)
                {
                    if (massFlux[f] > 0.)
                    {
                        gbc.applyExtrapolationBC(f);
                    }
                    else
                    {
                        gbc.applyDirichletBC(f,bVelocity[f]);
                    }
                }
                if (bc.bcType == "PressureBoundary")
                {
                    fixedPressureMomentumBC(faces,mesh,
                                            ls.getMatrix(), ls.getX(), ls.getB());
                }
            }
            else if ((bc.bcType == "Symmetry"))
            {
                gbc.applySymmetryBC();
            }
	   
            else
              throw CException(bc.bcType + " not implemented for FlowModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<VectorT3,DiagTensorT3,T> gbc(faces,mesh,
                                                    _geomFields,
                                                    _flowFields.velocity,
                                                    _flowFields.momentumFlux,
                                                    ls.getMatrix(), ls.getX(), ls.getB());

            gbc.applyInterfaceBC();
        }

    }
    DiscrList discretizations2;
    shared_ptr<Discretization>
      ud(new Underrelaxer<VectorT3,DiagTensorT3,T>
         (_meshes,_flowFields.velocity,
          _options["momentumURF"]));

    discretizations2.push_back(ud);

    linearizer.linearize(discretizations2,_meshes,ls.getMatrix(),
                         ls.getX(), ls.getB());

  }


  MFRPtr solveMomentum()
  {
    LinearSystem ls;

    initMomentumLinearization(ls);
    ls.initAssembly();
    linearizeMomentum(ls);
    ls.initSolve();
    

    // save current velocity for use in continuity discretization
    _previousVelocity = dynamic_pointer_cast<Field>(_flowFields.velocity.newCopy());
    
    //AMG solver(ls);
    MFRPtr rNorm = _options.getMomentumLinearSolver().solve(ls);
    
    if (!_initialMomentumNorm) _initialMomentumNorm = rNorm;
        
    _options.getMomentumLinearSolver().cleanup();
    
    ls.postSolve();
    
    ls.updateSolution();
   
    // save the momentum ap coeffficients for use in continuity discretization
    _momApField = shared_ptr<Field>(new Field("momAp"));
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
        const VVMatrix& vvMatrix =
          dynamic_cast<const VVMatrix&>(ls.getMatrix().getMatrix(vIndex,vIndex));
        const VVDiagArray& momAp = vvMatrix.getDiag();
        _momApField->addArray(cells,dynamic_pointer_cast<ArrayBase>(momAp.newCopy()));
    }
    _momApField->syncLocal();
    return rNorm;
  }



  void interfaceContinuityBC(const Mesh& mesh,
                             const StorageSite& faces,
                             MultiFieldMatrix& matrix,
                             MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
    
    const int nFaces = faces.getCount();

    for(int f=0; f<nFaces; f++)
    {
        int c0 = faceCells(f,0);
        int c1 = faceCells(f,1);
        // either c0 or c1 could be the exterior cell
        if (c1 >= cells.getSelfCount())
        {
            rCell[c1] = 0;
        }
        else
        {
            rCell[c0] = 0;
        }
    }
  }


  void correctVelocityBoundary(const Mesh& mesh,
                               const StorageSite& faces,
                               const MultiField& ppField) 
  {
    const StorageSite& cells = mesh.getCells();

    const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
    
    VectorT3Array& V = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
    const TArray& pp = dynamic_cast<const TArray&>(ppField[pIndex]);
    
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        const T ppFace = pp[c1];
        const VectorT3 ppA = ppFace*faceArea[f];

        V[c0] += ppA/momAp[c0];
    }
  }


  void correctMassFluxBoundary(const StorageSite& faces,
                               const MultiField& ppField)
  {
    MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);
    const TArray& dMassFlux = dynamic_cast<const TArray&>(ppField[mfIndex]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
    
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        massFlux[f] -= dMassFlux[f];
    }
  }

  void correctPressure(const Mesh& mesh,
                       const MultiField& xField)
  {
    const StorageSite& cells = mesh.getCells();
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    
    TArray& p = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
    const TArray& pp = dynamic_cast<const TArray&>(xField[pIndex]);

    const T pressureURF(_options["pressureURF"]);
      
    const int nCells = cells.getCountLevel1();
    for(int c=0; c<nCells; c++)
    {
        p[c] += pressureURF*(pp[c]-_referencePP);
    }
  }


  void correctVelocityExplicit(const Mesh& mesh,
                               const MultiField& xField)
  {
    const StorageSite& cells = mesh.getCells();
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    
    VectorT3Array& V = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
    const VectorT3Array& Vp = dynamic_cast<const VectorT3Array&>(xField[vIndex]);

    const T velocityURF(_options["velocityURF"]);
      
    const int nCells = cells.getCountLevel1();
    for(int c=0; c<nCells; c++)
    {
        V[c] += velocityURF*Vp[c];
    }

    // boundary
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        const StorageSite& faces = fg.site;
        MultiField::ArrayIndex fluxIndex(&_flowFields.momentumFlux,&faces);
        VectorT3Array& momFlux =
          dynamic_cast<VectorT3Array&>(_flowFields.momentumFlux[faces]);
        const VectorT3Array& dmomFlux =
          dynamic_cast<const VectorT3Array&>(xField[fluxIndex]);
        
        const int nFaces = faces.getCount();
        for(int f=0; f<nFaces; f++)
        {
            momFlux[f] += dmomFlux[f];
        }
    }
  }


  // set the first cell of the first mesh to be a Dirichlet point
  void setDirichlet(MultiFieldMatrix& mfmatrix,
                    MultiField& rField)
  {
    const Mesh& mesh = *_meshes[0];
    
    const StorageSite& cells = mesh.getCells();
    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
    
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(mfmatrix.getMatrix(pIndex,pIndex));

    PPDiagArray& ppDiag = ppMatrix.getDiag();
    PPDiagArray& ppCoeff = ppMatrix.getOffDiag();

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);

    const CRConnectivity& cr = ppMatrix.getConnectivity();

    const Array<int>& row = cr.getRow();

    const int nCells = cells.getSelfCount();

#ifdef FVM_PARALLEL
    const Array<int>& localToGlobal = mesh.getLocalToGlobal();
#endif
    
    _globalRefCellID = numeric_limits<int>::max();
    for ( int n = 0; n < nCells; n++ )
    {
        if ( ibType[n] == Mesh::IBTYPE_FLUID )
        {
#ifdef FVM_PARALLEL
            int glblIndx = localToGlobal[n];
#else
            int glblIndx = n;
#endif
              
            if ( glblIndx < _globalRefCellID )
	      _globalRefCellID = glblIndx;
        }
    } 
    
#ifdef FVM_PARALLEL
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &_globalRefCellID, 1, MPI::INT, MPI::MIN);

    
    _globalRefProcID = -1;
    const map<int,int>& globalToLocal = mesh.getGlobalToLocal();

    if ( globalToLocal.count(_globalRefCellID) > 0 ){
       _globalRefProcID = MPI::COMM_WORLD.Get_rank();
    }
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &_globalRefProcID, 1, MPI::INT, MPI::MAX);
    
    
    if ( globalToLocal.count(_globalRefCellID) > 0)
    {
       int nc =  globalToLocal.find(_globalRefCellID)->second;

#else
       int nc = _globalRefCellID;
#endif
       
       ppDiag[nc] = -1.;
       rCell[nc]=0.;
       for(int nb=row[nc]; nb<row[nc+1]; nb++)
         ppCoeff[nb] = 0;
#ifdef PV_COUPLED
       if (mfmatrix.hasMatrix(pIndex,vIndex))
       {
          PVMatrix& pvMatrix =
            dynamic_cast<PVMatrix&>(mfmatrix.getMatrix(pIndex,vIndex));
        
          PVDiagArray& pvDiag = pvMatrix.getDiag();
          PVDiagArray& pvCoeff = pvMatrix.getOffDiag();

          pvDiag[nc] = NumTypeTraits<VectorT3>::getZero();
          for(int nb=row[nc]; nb<row[nc+1]; nb++)
             pvCoeff[nb] = NumTypeTraits<VectorT3>::getZero();
       }
#endif
       
#ifdef FVM_PARALLEL
    }
#endif
  }

  

  void linearizeContinuity(LinearSystem& ls)
  {
    MultiFieldMatrix& matrix = ls.getMatrix();
    MultiField& x = ls.getX();
    MultiField& b = ls.getB();

    
    this->_useReferencePressure = true;
    const int numMeshes = _meshes.size();

    T netFlux(0.);
    

    _flowFields.pressureGradient.syncLocal();

    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const StorageSite& iFaces = mesh.getInteriorFaceGroup().site;
        
        MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
        
        // interior
        
        netFlux += discretizeMassFluxInterior(mesh,iFaces,matrix,x,b);
        
        // interfaces
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            netFlux += discretizeMassFluxInterior(mesh,faces,matrix,x,b);

            // this just sets the residual for the ghost cell to zero
            interfaceContinuityBC(mesh,faces,matrix,b);
        }
        
        // boundaries
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            const FlowBC<T>& bc = *_bcMap[fg.id];
            
            if ((bc.bcType == "NoSlipWall") ||
                (bc.bcType == "SlipJump") ||
                (bc.bcType == "Symmetry") ||
                (bc.bcType == "VelocityBoundary"))
            {
                
                netFlux += fixedFluxContinuityBC(faces,mesh,matrix,x,b,bc);
            }
            else if (bc.bcType == "PressureBoundary")
            {
                netFlux += fixedPressureContinuityBC(faces,mesh,matrix,x,b);
                this->_useReferencePressure=false;
            }	
            else		
              throw CException(bc.bcType + " not implemented for FlowModel");
            
            MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);
            typedef DiagonalMatrix<T,T> FFMatrix;
            FFMatrix& dFluxdFlux =
              dynamic_cast<FFMatrix&>(matrix.getMatrix(mfIndex,mfIndex));
            dFluxdFlux.unitize();
        }

	// add additional imbalance to netFlux due to transient term
	// when we have deforming mesh
        if (_geomFields.volumeN1.hasArray(cells))
        {

            // this block computes the net sum of d(rho)/dt over all the cells
            const TArray& density =
              dynamic_cast<const TArray&>(_flowFields.density[cells]);	
            const TArray& cellVolume =
              dynamic_cast<const TArray&>(_geomFields.volume[cells]);
            const int nCells = cells.getSelfCount();
            T _dT(_options["timeStep"]);
            T onePointFive(1.5);
            T two(2.0);
            T pointFive(0.5);

            const TArray& cellVolumeN1 = 
              dynamic_cast<const TArray&>(_geomFields.volumeN1[cells]);

            if (_geomFields.volumeN2.hasArray(cells))
	    {
                // second order time discretization

                const TArray& cellVolumeN2 = 
		  dynamic_cast<const TArray&>(_geomFields.volumeN2[cells]);
		for(int c=0; c<nCells; c++)
		{
		    const T rhobydT = density[c]/_dT;
		    const T term1 = onePointFive*cellVolume[c];
		    const T term2 = two*cellVolumeN1[c];
		    const T term3 = pointFive*cellVolumeN2[c];
		    netFlux -= rhobydT*(term1 - term2 + term3);
		}
	    }
            else
            {
                // first order time discretization
		for(int c=0; c<nCells; c++)
		{       
		    const T rhobydT = density[c]/_dT;
		    netFlux -= rhobydT*(cellVolume[c] - cellVolumeN1[c]);
		}
	    }

            // compute the net integral of the grid flux term 
            foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
            {
                const FaceGroup& fg = *fgPtr;
                const StorageSite& faces = fg.site;
                const int nFaces = faces.getCount();
                const CRConnectivity& faceCells = mesh.getFaceCells(faces);
                const VectorT3Array& faceArea =
                  dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
                const VectorT3Array& faceVel =
                  dynamic_cast<const VectorT3Array&>(_geomFields.faceVel[faces]);
                for(int f=0; f<nFaces; f++)
                {
                    const int c0 = faceCells(f,0);
                    netFlux += density[c0]*dot(faceVel[f],faceArea[f]);
                }
            }
        }
    }
  
 
      //sum netflux globalally MPI::
#ifdef FVM_PARALLEL
      int count = 1;
      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE, &netFlux, count, MPI::DOUBLE, MPI::SUM);

      int useReferencePressureInt = int(this->_useReferencePressure );
      MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &useReferencePressureInt, count, MPI::INT, MPI::PROD);
      this->_useReferencePressure = bool(useReferencePressureInt);
#endif
        //cout << "net boundary flux = " << netFlux << endl;
    
    if (this->_useReferencePressure)
    {
        // we have all fixed mass flux boundary problem. In case the
        // net mass flux sum over all the boundaries does not add up
        // to zero we distribute that over the cells as a volumetric
        // as a matching source/sink so that the continuity equation
        // can converge to a zero residual. In case of immersed
        // boundary problems we need to use only the fluid cells since
        // all the other cells will always have zero corrections.
        
        T volumeSum(0.);
        
        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            const StorageSite& cells = mesh.getCells();
            const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
            for(int c=0; c<cells.getSelfCount(); c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
                volumeSum += cellVolume[c];
        }

#ifdef FVM_PARALLEL
        MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &volumeSum, 1, MPI::DOUBLE, MPI::SUM);
#endif	

        netFlux /= volumeSum;

        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            const StorageSite& cells = mesh.getCells();
            const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
            
            MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
            TArray& rCell = dynamic_cast<TArray&>(b[pIndex]);
            
            for(int c=0; c<cells.getSelfCount(); c++)
            {
              if (ibType[c] == Mesh::IBTYPE_FLUID)
                rCell[c] += netFlux*cellVolume[c];
            }
        }

        // when all the mass fluxes are specified the continuity
        // equation also has an extra degree of freedom. This sets the
        // first cell to have a zero correction to account for this.
        

      setDirichlet(matrix,b);

    }
  }

  void setReferencePP(const MultiField& ppField)
  {
    if (_useReferencePressure)
    {
        const Mesh& mesh = *_meshes[0];
        
        const StorageSite& cells = mesh.getCells();
        
        MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
        const TArray& pp = dynamic_cast<const TArray&>(ppField[pIndex]);

#ifndef FVM_PARALLEL
        _referencePP = pp[_globalRefCellID];  
#endif

        
#ifdef FVM_PARALLEL    
        _referencePP = 0.0;
        const map<int,int>& globalToLocal = mesh.getGlobalToLocal();
	if ( MPI::COMM_WORLD.Get_rank() == _globalRefProcID ){
	    int localID = globalToLocal.find(_globalRefCellID)->second;
           _referencePP = pp[localID];
        }	   
	   
        int count = 1;
	MPI::COMM_WORLD.Bcast( &_referencePP, count, MPI::DOUBLE, _globalRefProcID);			

#endif
       //broadcast this value (referencePP) to all 

    }
    else
      _referencePP = 0.0;
  }

  void computeContinuityResidual()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        const StorageSite& faces = mesh.getFaces();

        TArray& r = dynamic_cast<TArray&>(_flowFields.continuityResidual[cells]);
        const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);

        const CRConnectivity& faceCells = mesh.getAllFaceCells();
        const int nFaces = faces.getCount();

        r.zero();
        for(int f=0; f<nFaces; f++)
        {
	  const int c0 = faceCells(f,0);
	  const int c1 = faceCells(f,1);

            r[c0] += massFlux[f];
            r[c1] -= massFlux[f];
        }
    }
  }

  void postContinuitySolve(LinearSystem& ls)
  {
    MultiFieldMatrix& matrix = ls.getMatrix();
    MultiField& ppField = ls.getDelta();

    setReferencePP(ppField);
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        const StorageSite& iFaces = mesh.getInteriorFaceGroup().site;
        
        MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
        
        MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
        const bool coupled =  matrix.hasMatrix(pIndex,vIndex);
        
        correctPressure(mesh,ppField);
        
        // interior
        
        correctMassFluxInterior(mesh,iFaces,matrix,ppField);
        
        if (coupled)
          correctVelocityExplicit(mesh,ppField);
        else if (_options.correctVelocity)
          correctVelocityInterior(mesh,iFaces,ppField);
        
        updateFacePressureInterior(mesh,iFaces);

        // interfaces
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            correctMassFluxInterior(mesh,faces,matrix,ppField);
            
            if (coupled)
              correctVelocityExplicit(mesh,ppField);
            else if (_options.correctVelocity)
              correctVelocityInterior(mesh,faces,ppField);
            
            updateFacePressureInterior(mesh,faces);
        }
        
        // boundary
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            const FlowBC<T>& bc = *_bcMap[fg.id];
              
            correctMassFluxBoundary(faces,ppField);
            
            if (!coupled && _options.correctVelocity)
              correctVelocityBoundary(mesh,faces,ppField);

            if (bc.bcType == "PressureBoundary")
            {
                T bp = bc["specifiedPressure"];
                pressureBoundaryPostContinuitySolve(faces,mesh,bp);
            }
            
            updateFacePressureBoundary(mesh,faces);
        }

    }

    _flowFields.velocity.syncLocal();
    _flowFields.pressure.syncLocal();

    computeContinuityResidual();
    
  }
  
  void initContinuityLinearization(LinearSystem& ls)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);

        ls.getX().addArray(pIndex,_flowFields.pressure.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();
        
        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        ls.getMatrix().addMatrix(pIndex,pIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_flowFields.massFlux,&faces);
            ls.getX().addArray(fIndex,_flowFields.massFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,pIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_flowFields.massFlux,&faces);
            ls.getX().addArray(fIndex,_flowFields.massFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            ls.getMatrix().addMatrix(fIndex,pIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            ls.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
    }
  }

  shared_ptr<LinearSystem> discretizeContinuity()
  {
    shared_ptr<LinearSystem> ls(new LinearSystem());

    initContinuityLinearization(*ls);
        
    ls->initAssembly();
    
    linearizeContinuity(*ls);

    ls->initSolve();

    return ls;
  }

  
  MFRPtr solveContinuity()
  {
    shared_ptr<LinearSystem> ls(discretizeContinuity());

    // discard previous velocity
    _previousVelocity = shared_ptr<Field>();

    MFRPtr rNorm = _options.getPressureLinearSolver().solve(*ls);

    if (!_initialContinuityNorm) _initialContinuityNorm = rNorm;

    ls->postSolve();
    _options.pressureLinearSolver->cleanup();
    
    postContinuitySolve(*ls);

    // discard the momentum ap coeffficients
    
    _momApField = shared_ptr<Field>();
    return rNorm;
  }


  bool advance(const int niter)
  {

    for(int n=0; n<niter; n++)
    { 
        MFRPtr mNorm = solveMomentum();
        MFRPtr cNorm = solveContinuity();

        if (_niters < 5)
        {
            _initialMomentumNorm->setMax(*mNorm);
            _initialContinuityNorm->setMax(*cNorm);
        }
        
        MFRPtr mNormRatio((*mNorm)/(*_initialMomentumNorm));
	MFRPtr cNormRatio((*cNorm)/(*_initialContinuityNorm));
        
#ifdef FVM_PARALLEL
	if ( MPI::COMM_WORLD.Get_rank() == 0 ){ // only root process
        if (_options.printNormalizedResiduals)
          {cout << _niters << ": " << *mNormRatio << ";" << *cNormRatio <<  endl;}
        else{
          cout << _niters << ": " << *mNorm << ";" << *cNorm <<  endl;}}
#endif

#ifndef FVM_PARALLEL
        if (_options.printNormalizedResiduals)
          {cout << _niters << ": " << *mNormRatio << ";" << *cNormRatio <<  endl;}
        else
          {cout << _niters << ": " << *mNorm << ";" << *cNorm <<  endl;}
#endif

        _niters++;
        if ((*mNormRatio < _options.momentumTolerance) &&
            (*cNormRatio < _options.continuityTolerance))
          return true;
    }
    return false;
  }

#ifdef PV_COUPLED
  bool advanceCoupled(const int niter)
  {
    const int numMeshes = _meshes.size();
    for(int n=0; n<niter; n++)
    { 
        LinearSystem ls;

        ls.setCoarseningField(_flowFields.pressure);
        initMomentumLinearization(ls);
        initContinuityLinearization(ls);
        
        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            
            const StorageSite& cells = mesh.getCells();
            MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
            MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
            
            const CRConnectivity& cellCells = mesh.getCellCells();
            
            shared_ptr<Matrix> mvp(new VPMatrix(cellCells));
            ls.getMatrix().addMatrix(vIndex,pIndex,mvp);

            shared_ptr<Matrix> mpv(new PVMatrix(cellCells));
            ls.getMatrix().addMatrix(pIndex,vIndex,mpv);
        }

        ls.initAssembly();

        linearizeMomentum(ls);

        // save current velocity for use in continuity discretization
        _previousVelocity = dynamic_pointer_cast<Field>(_flowFields.velocity.newCopy());
        
        // save the momentum ap coeffficients for use in continuity discretization
        _momApField = shared_ptr<Field>(new Field("momAp"));
        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            
            const StorageSite& cells = mesh.getCells();
            MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
            const VVMatrix& vvMatrix =
              dynamic_cast<const VVMatrix&>(ls.getMatrix().getMatrix(vIndex,vIndex));
            const VVDiagArray& momAp = vvMatrix.getDiag();
            _momApField->addArray(cells,dynamic_pointer_cast<ArrayBase>(momAp.newCopy()));
        }
        
        linearizeContinuity(ls);
        
        ls.initSolve();

        MFRPtr rNorm = _options.coupledLinearSolver->solve(ls);

        if (!_initialCoupledNorm) _initialCoupledNorm = rNorm;
        
        ls.postSolve();

        postContinuitySolve(ls);
        
        _options.coupledLinearSolver->cleanup();
    
        if (_niters < 5)
          _initialCoupledNorm->setMax(*rNorm);
        
        MFRPtr normRatio((*rNorm)/(*_initialCoupledNorm));
        
        if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *normRatio <<  endl;
        else
          cout << _niters << ": " << *rNorm <<  endl;

        _niters++;

        _momApField = shared_ptr<Field>();

        if (*normRatio < _options.momentumTolerance)
          return true;
    }
    return false;
  }
#endif
  
#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
  
  void dumpContinuityMatrix(const string fileBase)
  {
    solveMomentum();
    shared_ptr<LinearSystem> ls(discretizeContinuity());

    MultiFieldMatrix& matrix = ls->getMatrix();
    MultiField& b = ls->getB();

    const Mesh& mesh = *_meshes[0];
    const StorageSite& cells = mesh.getCells();
    
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(matrix.getMatrix(pIndex,pIndex));

    PPDiagArray& ppDiag = ppMatrix.getDiag();
    PPDiagArray& ppCoeff = ppMatrix.getOffDiag();

    TArray& rCell = dynamic_cast<TArray&>(b[pIndex]);

    const CRConnectivity& cr = ppMatrix.getConnectivity();

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
    
    string matFileName = fileBase + ".mat";
    FILE *matFile = fopen(matFileName.c_str(),"wb");
    
    fprintf(matFile,"%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matFile,"%d %d %d\n", nCells,nCells,nCoeffs);

    for(int i=0; i<nCells; i++)
    {
        fprintf(matFile,"%d %d %lf\n", i+1, i+1, ppDiag[i]);
        for(int jp=row[i]; jp<row[i+1]; jp++)
        {
            const int j = col[jp];
            if (j<nCells)
              fprintf(matFile,"%d %d %lf\n", i+1, j+1, ppCoeff[jp]);
        }
    }

    fclose(matFile);

    string rhsFileName = fileBase + ".rhs";
    FILE *rhsFile = fopen(rhsFileName.c_str(),"wb");
    
    for(int i=0; i<nCells; i++)
      fprintf(rhsFile,"%lf\n",-rCell[i]);

    fclose(rhsFile);
  }
#endif
  
  void printBCs()
  {
    foreach(typename FlowBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename FlowBC<T>::value_type& vp, *pos.second)
        {
            cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
        }
    }
  }

  VectorT3 getPressureIntegral(const Mesh& mesh, const int faceGroupId)
  {
    VectorT3 r(VectorT3::getZero());
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
            const int nFaces = faces.getCount();
            const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
            for(int f=0; f<nFaces; f++)
              r += faceArea[f]*facePressure[f];

            found=true;
        }
    }
    if (!found)
      throw CException("getPressureIntegral: invalid faceGroupID");
    return r;
  }

  VectorT3 getPressureIntegralonIBFaces(const Mesh& mesh)
  {
    VectorT3 r(VectorT3::getZero());
    const StorageSite& ibFaces = mesh.getIBFaces();
    const StorageSite& faces = mesh.getFaces();
    const StorageSite& cells = mesh.getCells();
    const Array<int>& ibFaceIndices = mesh.getIBFaceList();
    const int nibf = ibFaces.getCount();
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
    const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);

    for ( int f = 0; f < nibf; f ++)
    {
        const int ibFaceIndex = ibFaceIndices[f];
        const int c0 = faceCells(ibFaceIndex,0);

        // need to check whether the face points in or out of the
        // fluid cell and adjust sign of area accordingly
        
        if (ibType[c0] == Mesh::IBTYPE_FLUID)
          r += faceArea[ibFaceIndex]*facePressure[ibFaceIndex];
        else
          r -= faceArea[ibFaceIndex]*facePressure[ibFaceIndex];
    }

#if FVM_PARALLEL    
     MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, r.getData(), 3, MPI::DOUBLE, MPI::SUM);
#endif     


    return r;
  }

 
  VectorT3 getMomentumFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    VectorT3 r(VectorT3::getZero());
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const VectorT3Array& momFlux =
              dynamic_cast<const VectorT3Array&>(_flowFields.momentumFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += momFlux[f];
            found=true;
        }
    }
    if (!found)
      throw CException("getMomentumFluxIntegral: invalid faceGroupID");
    return r;
  }

  VectorT3 getMomentumDerivativeIntegral(const Mesh& mesh)
  {
    VectorT3 r(VectorT3::getZero());
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();

    const TArray& density =
          dynamic_cast<const TArray&>(_flowFields.density[cells]);
    const VectorT3Array& v =
          dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
    const VectorT3Array& vN1 =
      dynamic_cast<const VectorT3Array&>(_flowFields.velocityN1[cells]);

    const TArray& volume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const T deltaT = _options["timeStep"];
    
    if (_flowFields.velocityN2.hasArray(cells))
    {
        // second order
        const VectorT3Array& vN2 =
          dynamic_cast<const VectorT3Array&>(_flowFields.velocityN2[cells]);
        T onePointFive(1.5);
        T two(2.0);
        T pointFive(0.5);
        
        for(int c=0; c<nCells; c++)
        {
            const T rhoVbydT = density[c]*volume[c]/deltaT;
            r += rhoVbydT*(onePointFive*v[c]- two*vN1[c]
                           + pointFive*vN2[c]);
        }
    }
    else
    {
        for(int c=0; c<nCells; c++)
        {
            const T rhoVbydT = density[c]*volume[c]/deltaT;
            r += rhoVbydT*(v[c]- vN1[c]);
        }
    }
    return r;
  }

  boost::shared_ptr<ArrayBase> getStressTensor(const Mesh& mesh, const ArrayBase& gcellIds)
  {
    typedef Array<StressTensor<T> > StressTensorArray;
    
    const StorageSite& cells = mesh.getCells();
    
    const Array<int>& cellIds = dynamic_cast<const Array<int> &>(gcellIds);
    const int nCells = cellIds.getLength();

    _velocityGradientModel.compute();

    const VGradArray& vGrad =
      dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);

    const TArray& pCell =
      dynamic_cast<const TArray&>(_flowFields.pressure[cells]);

    const TArray& mu = dynamic_cast<const TArray&>(getViscosityField()[cells]);

    boost::shared_ptr<StressTensorArray> stressTensorPtr( new StressTensorArray(nCells));
    StressTensorArray& stressTensor = *stressTensorPtr;

    for(int n=0; n<nCells; n++)
    {
        const int c = cellIds[n];
        const VGradType& vg = vGrad[c];
        VGradType vgPlusTranspose = vGrad[c];

        for(int i=0;i<3;i++)
          for(int j=0;j<3;j++)
            vgPlusTranspose[i][j] += vg[j][i];
        
        stressTensor[n][0] = vgPlusTranspose[0][0]*mu[c] - pCell[c];
        stressTensor[n][1] = vgPlusTranspose[1][1]*mu[c] - pCell[c];
        stressTensor[n][2] = vgPlusTranspose[2][2]*mu[c] - pCell[c];
        stressTensor[n][3] = vgPlusTranspose[0][1]*mu[c];
        stressTensor[n][4] = vgPlusTranspose[1][2]*mu[c];
        stressTensor[n][5] = vgPlusTranspose[2][0]*mu[c];
    }

    return stressTensorPtr;
  }

  void ComputeStressTensorES(const StorageSite& solidFaces)
  {
    typedef Array<StressTensor<T> > StressTensorArray;
    const int nSolidFaces = solidFaces.getCount();


    //interface
    VectorT3Array& IntVel =
      dynamic_cast<VectorT3Array&>(_flowFields.InterfaceVelocity[solidFaces]);
    TArray& IntPress =
      dynamic_cast<TArray&>(_flowFields.InterfacePressure[solidFaces]); 
    TArray& IntDens =
      dynamic_cast<TArray&>(_flowFields.InterfaceDensity[solidFaces]);
   
    //const TArray& mu = dynamic_cast<const TArray&>(getViscosityField()[cells]);   
    StressTensorArray& stressTensor = dynamic_cast<StressTensorArray &>(_flowFields.InterfaceStress[solidFaces]);
   

    const int numMeshes = _meshes.size();
 
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	//const int nCells = cells.getCount();
	
	_velocityGradientModel.compute();

	const VGradArray& vGrad =
	  dynamic_cast<VGradArray&>(_flowFields.velocityGradient[cells]);
	const TArray& pCell =
	  dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
	const TArray& dCell =
	  dynamic_cast<const TArray&>(_flowFields.density[cells]);
    	const VectorT3Array& vCell =
	  dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);


	const CRConnectivity& _faceCells=mesh.getFaceCells(solidFaces);
	
	for(int f=0; f<nSolidFaces; f++)
	  {
	   
	    const int c0 = _faceCells(f,0);  
	    const int c1 = _faceCells(f,1);     
	    
	    const VGradType& vg = vGrad[c0];

	    VGradType vgPlusTranspose = vGrad[c0];
	    
	    for(int i=0;i<3;i++)
	      for(int j=0;j<3;j++)
		vgPlusTranspose[i][j] += vg[j][i];
	    
	    stressTensor[f][0] = vgPlusTranspose[0][0];
	    stressTensor[f][1] = vgPlusTranspose[1][1];
	    stressTensor[f][3] = vgPlusTranspose[0][1];
	    stressTensor[f][4] = vgPlusTranspose[1][2];
	    stressTensor[f][5] = vgPlusTranspose[2][0];
	    
	    //copy density pressure and velocity into Interface arrays
	    IntDens[f]=dCell[c1];
	    IntPress[f]=pCell[c1];
	    IntVel[f][0]=vCell[c1][0];
	    IntVel[f][1]=vCell[c1][1];
	    IntVel[f][2]=vCell[c1][2];
	      }
      }
    
    
  }
  

  void getTraction(const Mesh& mesh)
  {
    const StorageSite& cells = mesh.getCells();

    const int nCells = cells.getSelfCount();

    shared_ptr<VectorT3Array> tractionXPtr(new VectorT3Array(nCells));
    tractionXPtr->zero();
    _flowFields.tractionX.addArray(cells,tractionXPtr);
    VectorT3Array& tractionX = *tractionXPtr;

    shared_ptr<VectorT3Array> tractionYPtr(new VectorT3Array(nCells));
    tractionYPtr->zero();
    _flowFields.tractionY.addArray(cells,tractionYPtr);
    VectorT3Array& tractionY = *tractionYPtr;

    shared_ptr<VectorT3Array> tractionZPtr(new VectorT3Array(nCells));
    tractionZPtr->zero();
    _flowFields.tractionZ.addArray(cells,tractionZPtr);
    VectorT3Array& tractionZ = *tractionZPtr;

    _velocityGradientModel.compute();

    const VGradArray& vGrad =
      dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);

    const TArray& pCell =
      dynamic_cast<const TArray&>(_flowFields.pressure[cells]);

    const TArray& mu = dynamic_cast<const TArray&>(getViscosityField()[cells]);
      
    for(int n=0; n<nCells; n++)
    {
        const VGradType& vg = vGrad[n];
        VGradType vgPlusTranspose = vGrad[n];

        for(int i=0;i<3;i++)
          for(int j=0;j<3;j++)
            vgPlusTranspose[i][j] += vg[j][i];
          
        tractionX[n][0] = vgPlusTranspose[0][0]*mu[n] - pCell[n];
        tractionX[n][1] = vgPlusTranspose[0][1]*mu[n];
        tractionX[n][2] = vgPlusTranspose[0][2]*mu[n];

        tractionY[n][0] = vgPlusTranspose[1][0]*mu[n];
        tractionY[n][1] = vgPlusTranspose[1][1]*mu[n] - pCell[n];
        tractionY[n][2] = vgPlusTranspose[1][2]*mu[n];

        tractionZ[n][0] = vgPlusTranspose[2][0]*mu[n];
        tractionZ[n][1] = vgPlusTranspose[2][1]*mu[n];
        tractionZ[n][2] = vgPlusTranspose[2][2]*mu[n] - pCell[n];
    }
  }

  void
  computeSolidSurfaceForce(const StorageSite& solidFaces, bool perUnitArea)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;

    const int nSolidFaces = solidFaces.getCount();

    _velocityGradientModel.compute();

    boost::shared_ptr<VectorT3Array>
      forcePtr( new VectorT3Array(nSolidFaces));
    VectorT3Array& force = *forcePtr;

    force.zero();
    _flowFields.force.addArray(solidFaces,forcePtr);

    const VectorT3Array& solidFaceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[solidFaces]);

    const TArray& solidFaceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[solidFaces]);
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const VGradArray& vGrad =
          dynamic_cast<const VGradArray&>(_flowFields.velocityGradient[cells]);

        const TArray& pCell =
          dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
        
        const TArray& mu =
          dynamic_cast<const TArray&>(getViscosityField()[cells]);
        
        //const FlowVC<T>& vc = *_vcMap[mesh.getID()];
            
        const CRConnectivity& solidFacesToCells
          = mesh.getConnectivity(solidFaces,cells);
        
        const IntArray& sFCRow = solidFacesToCells.getRow();
        const IntArray& sFCCol = solidFacesToCells.getCol();

        GeomFields::SSPair key1(&solidFaces,&cells);
        const IMatrix& mIC =
          dynamic_cast<const IMatrix&>
          (*_geomFields._interpolationMatrices[key1]);

        const Array<T>& iCoeffs = mIC.getCoeff();
        
        for(int f=0; f<nSolidFaces; f++)
        {
            StressTensor<T> stress = NumTypeTraits<StressTensor<T> >::getZero();
            for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
            {
                const int c = sFCCol[nc];
                const VGradType& vg = vGrad[c];
                VGradType vgPlusTranspose = vGrad[c];
                
                for(int i=0;i<3;i++)
                  for(int j=0;j<3;j++)
                    vgPlusTranspose[i][j] += vg[j][i];

                const T coeff = iCoeffs[nc];
               
                stress[0] += coeff*(vgPlusTranspose[0][0]*mu[c] - pCell[c]);
                stress[1] += coeff*(vgPlusTranspose[1][1]*mu[c] - pCell[c]);
                stress[2] += coeff*(vgPlusTranspose[2][2]*mu[c] - pCell[c]);
                stress[3] += coeff*(vgPlusTranspose[0][1]*mu[c]);
                stress[4] += coeff*(vgPlusTranspose[1][2]*mu[c]);
                stress[5] += coeff*(vgPlusTranspose[2][0]*mu[c]);
            }

            const VectorT3& Af = solidFaceArea[f];
            force[f][0] = Af[0]*stress[0] + Af[1]*stress[3] + Af[2]*stress[5];
            force[f][1] = Af[0]*stress[3] + Af[1]*stress[1] + Af[2]*stress[4];
            force[f][2] = Af[0]*stress[5] + Af[1]*stress[4] + Af[2]*stress[2];
            if (perUnitArea)
            {
                force[f] /= solidFaceAreaMag[f];
            }
        }
    }
  }


  VectorT3 getMomentumFluxIntegralonIBFaces(const Mesh& mesh)
  {
    VectorT3 r(VectorT3::getZero());
    const StorageSite& ibFaces = mesh.getIBFaces();
    const StorageSite& faces = mesh.getFaces();
    const Array<int>& ibFaceIndices = mesh.getIBFaceList();
    const int nibf = ibFaces.getCount();
    const VectorT3Array& momFlux =
              dynamic_cast<const VectorT3Array&>(_flowFields.momentumFlux[faces]);
    for ( int f = 0; f < nibf; f ++){
      const int ibFaceIndex = ibFaceIndices[f];
      r += momFlux[ibFaceIndex];
    }
    if (nibf == 0)
      throw CException("getMomentumFluxIntegralonIBFaces:  no ibFaces found!");
    return r;
  }
  

  void printPressureIntegrals()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            
            VectorT3 r(VectorT3::getZero());
            
            const StorageSite& faces = fg.site;
            const VectorT3Array& faceArea =
              dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
            const int nFaces = faces.getCount();
            const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
            for(int f=0; f<nFaces; f++)
              r += faceArea[f]*facePressure[f];

            cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
        }
    }
  }
  
  void printMomentumFluxIntegrals()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            
            VectorT3 r(VectorT3::getZero());
            
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const VectorT3Array& momFlux =
              dynamic_cast<const VectorT3Array&>(_flowFields.momentumFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += momFlux[f];

            cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
        }
    }
  }
  
  void printMassFluxIntegrals()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            
            T r(0.);
            
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();
            const TArray& massFlux = dynamic_cast<const TArray&>(_flowFields.massFlux[faces]);
            for(int f=0; f<nFaces; f++)
              r += massFlux[f];

            cout << "Mesh " << mesh.getID() << " faceGroup " << fg.id << " : " << r <<  endl;
        }
    }
  }
  
#include "FlowModelPressureBC.h"
#include "FlowModelVelocityBC.h"
#include "FlowModelInterior.h"
#include "FlowModelSlipJump.h"
  
private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  FlowFields& _flowFields;
  FlowBCMap _bcMap;
  FlowVCMap _vcMap;
  
  FlowModelOptions<T> _options;
  GradientModel<VectorT3> _velocityGradientModel;
  GradientModel<T> _pressureGradientModel;
  
  MFRPtr _initialMomentumNorm;
  MFRPtr _initialContinuityNorm;
  MFRPtr _initialCoupledNorm;
  int _niters;

  shared_ptr<Field> _previousVelocity;
  shared_ptr<Field> _momApField;

  bool _useReferencePressure;
  int  _globalRefCellID;
  int  _globalRefProcID;
  T _referencePP;
  //AMG _momSolver;
  //AMG _continuitySolver;
  map<string,shared_ptr<ArrayBase> > _persistenceData;
};

template<class T>
FlowModel<T>::FlowModel(const GeomFields& geomFields,
                        FlowFields& thermalFields,
                        const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,thermalFields,meshes))
{
  logCtor();
}


template<class T>
FlowModel<T>::~FlowModel()
{
  logDtor();
}

template<class T>
void
FlowModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename FlowModel<T>::FlowBCMap&
FlowModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename FlowModel<T>::FlowVCMap&
FlowModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
FlowModelOptions<T>&
FlowModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
FlowModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
bool
FlowModel<T>::advance(const int niter)
{
  return _impl->advance(niter);
}

#ifdef PV_COUPLED
template<class T>
bool
FlowModel<T>::advanceCoupled(const int niter)
{
  return _impl->advanceCoupled(niter);
}
#endif

#if 0
template<class T>
LinearSolver&
FlowModel<T>::getMomentumSolver()
{
  return _impl->getMomentumSolver();
}

template<class T>
LinearSolver&
FlowModel<T>::getContinuitySolver()
{
  return _impl->getContinuitySolver();
}
#endif

template<class T>
void
FlowModel<T>::updateTime()
{
  _impl->updateTime();
}

template<class T>
void
FlowModel<T>::printPressureIntegrals()
{
  _impl->printPressureIntegrals();
}

template<class T>
void
FlowModel<T>::printMomentumFluxIntegrals()
{
  _impl->printMomentumFluxIntegrals();
}

template<class T>
void
FlowModel<T>::printMassFluxIntegrals()
{
  _impl->printMassFluxIntegrals();
}

template<class T>
Vector<T,3>
FlowModel<T>::getPressureIntegral(const Mesh& mesh, const int faceGroupId)
{
 return  _impl->getPressureIntegral(mesh,faceGroupId);
}

template<class T>
Vector<T,3>
FlowModel<T>::getMomentumFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
 return  _impl->getMomentumFluxIntegral(mesh,faceGroupId);
}

template<class T>
Vector<T,3>
FlowModel<T>::getMomentumDerivativeIntegral(const Mesh& mesh)
{
 return  _impl->getMomentumDerivativeIntegral(mesh);
}

template<class T>
Vector<T,3>
FlowModel<T>::getPressureIntegralonIBFaces(const Mesh& mesh)
{
 return  _impl->getPressureIntegralonIBFaces(mesh);
}

template<class T>
Vector<T,3>
FlowModel<T>::getMomentumFluxIntegralonIBFaces(const Mesh& mesh)
{
 return  _impl->getMomentumFluxIntegralonIBFaces(mesh);
}

template<class T>
void
FlowModel<T>::computeIBFaceVelocity(const StorageSite& particles)
{
  return _impl->computeIBFaceVelocity(particles);
}

template<class T>
void
FlowModel<T>::computeSolidSurfaceForce(const StorageSite& particles)
{
  return _impl->computeSolidSurfaceForce(particles,false);
}

template<class T>
void
FlowModel<T>::computeSolidSurfaceForcePerUnitArea(const StorageSite& particles)
{
  return _impl->computeSolidSurfaceForce(particles,true);
}


template<class T>
boost::shared_ptr<ArrayBase>
FlowModel<T>::getStressTensor(const Mesh& mesh, const ArrayBase& cellIds)
{
  return _impl->getStressTensor(mesh, cellIds);
}

template<class T>
void
FlowModel<T>::getTraction(const Mesh& mesh)
{
  return  _impl->getTraction(mesh);
}

template<class T>
void
FlowModel<T>:: ComputeStressTensorES(const StorageSite& particles)
{
  return  _impl-> ComputeStressTensorES(particles);
}
template<class T>
map<string,shared_ptr<ArrayBase> >&
FlowModel<T>::getPersistenceData()
{
  return _impl->getPersistenceData();
}

template<class T>
void
FlowModel<T>::restart()
{
  _impl->restart();
}
  



#if !(defined(USING_ATYPE_TANGENT) || defined(USING_ATYPE_PC))
template<class T>  
void
FlowModel<T>::dumpContinuityMatrix(const string fileBase)
{
  _impl->dumpContinuityMatrix(fileBase);
}

#endif
