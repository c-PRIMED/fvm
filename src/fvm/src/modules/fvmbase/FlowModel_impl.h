#include "Mesh.h"

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
#include "IbmDiscretization.h"

#include "Underrelaxer.h"
#include "MomentumPressureGradientDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"
#include "MomentumIBDiscretization.h"
#include "StressTensor.h"

template<class T>
class FlowModel<T>::Impl
{
public:
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

  typedef CRMatrixRect<VectorT3T,VectorT3,T> PVMatrix;
  typedef typename PVMatrix::DiagArray PVDiagArray;
  typedef typename PVMatrix::PairWiseAssembler PVAssembler;

  typedef CRMatrixRect<VectorT3,T,VectorT3> VPMatrix;
  typedef typename VPMatrix::DiagArray VPDiagArray;
  typedef typename VPMatrix::PairWiseAssembler VPAssembler;


  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;

  typedef Gradient<T> PGradType;
  typedef Array<PGradType> PGradArray;

  typedef FluxJacobianMatrix<T,T> FMatrix;
  
  Impl(const GeomFields& geomFields,
       FlowFields& thermalFields,
       const MeshList& meshes) :
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

        shared_ptr<VectorT3Array> vCell(new VectorT3Array(cells.getCount()));

        VectorT3 initialVelocity;
        initialVelocity[0] = _options["initialXVelocity"];
        initialVelocity[1] = _options["initialYVelocity"];
        initialVelocity[2] = _options["initialZVelocity"];
        *vCell = initialVelocity;

        _flowFields.velocity.addArray(cells,vCell);

        if (_options.transient)
        {
            _flowFields.velocityN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(vCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
              _flowFields.velocityN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(vCell->newCopy()));

        }
        
        shared_ptr<TArray> pCell(new TArray(cells.getCount()));
        shared_ptr<TArray> pFace(new TArray(faces.getCount()));
        *pCell = _options["initialPressure"];
        *pFace = _options["initialPressure"];
        _flowFields.pressure.addArray(cells,pCell);
        _flowFields.pressure.addArray(faces,pFace);


        shared_ptr<TArray> rhoCell(new TArray(cells.getCount()));
        *rhoCell = vc["density"];
        _flowFields.density.addArray(cells,rhoCell);

        shared_ptr<TArray> muCell(new TArray(cells.getCount()));
        *muCell = vc["viscosity"];
        _flowFields.viscosity.addArray(cells,muCell);

        shared_ptr<PGradArray> gradp(new PGradArray(cells.getCount()));
        gradp->zero();
        _flowFields.pressureGradient.addArray(cells,gradp);

        shared_ptr<TArray> ci(new TArray(cells.getCount()));
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
	// debug use
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
	  double r = mag(faceCentroid[fID]-center);
	  double angle = atan2(faceCentroid[fID][1]-center[1],faceCentroid[fID][0]-center[0]);
	  //(*ibV)[f][0]=-angV*r*sin(angle);
	  //(*ibV)[f][1]=angV*r*cos(angle);
	  //(*ibV)[f][2]=0.0;
	  (*ibV)[f][0]=0.001;
	  (*ibV)[f][1]=0.0;
	  (*ibV)[f][2]=0.0;
	}
	  
	//for(int f=0; f<ibFaces.getCount();f++){
	//cout<<f<<" "<<(*ibV)[f]<<endl;
	//}

#endif
        _flowFields.velocity.addArray(ibFaces,ibV);
    }

  }


 void computeIBandSolidVelocity(const StorageSite& particles)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    typedef CRMatrixTranspose<T,VectorT3,VectorT3> IMatrixV3;

    const VectorT3Array& pV =
      dynamic_cast<const VectorT3Array&>(_flowFields.velocity[particles]);
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
     	const int nCells = cells.getCount();
	
	GeomFields::SSPair key2(&cells,&particles);
        const IMatrix& mIP =
          dynamic_cast<const IMatrix&>
          (*_geomFields._interpolationMatrices[key2]);

        IMatrixV3 mIPV(mIP);


        shared_ptr<VectorT3Array> icV(new VectorT3Array(nCells));

        icV->zero();

       	mIPV.multiplyAndAdd(*icV,pV);
	
	//only modify the solid cell and IB cell velocity by interpolating of particles
	VectorT3Array& cV =
	  dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
		
	for (int c = 0; c < nCells; c ++){
	  const int cellType = mesh.getIBTypeForCell(c);
	  if (cellType == Mesh::IBTYPE_REALBOUNDARY){
	    cV[c] = (*icV)[c];
	  }
	}
	
	//_flowFields.velocity.addArray(cells,cV);
    }
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

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,VectorT3>(faceCells));
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

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,VectorT3>(faceCells));
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
          _flowFields.viscosity,
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

    shared_ptr<Discretization>
      ud(new Underrelaxer<VectorT3,DiagTensorT3,T>
         (_meshes,_flowFields.velocity,
          _options["momentumURF"]));

    discretizations.push_back(dd);
    discretizations.push_back(cd);
    discretizations.push_back(pd);
    discretizations.push_back(ud);

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

#if 0
    shared_ptr<Discretization>
      id(new IbmDiscretization<VectorT3,DiagTensorT3,T>
         (_meshes,_geomFields, _flowFields));

    discretizations.push_back(id);
#endif

    shared_ptr<Discretization>
      ibm(new MomentumIBDiscretization<VectorT3,DiagTensorT3,T>
             (_meshes,_geomFields,_flowFields));
      
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
                                                    ls.getMatrix(), ls.getX(), ls.getB());

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
            else if (bc.bcType == "VelocityBoundary")
                     //(bc.bcType == "PressureBoundary"))
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
            }
            else if (bc.bcType == "PressureBoundary")
            {
            }
            else if ((bc.bcType == "Symmetry"))
            {
                VectorT3 zeroFlux(NumTypeTraits<VectorT3>::getZero());
                gbc.applyNeumannBC(zeroFlux);
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
      
    const int nCells = cells.getCount();
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
      
    const int nCells = cells.getCount();
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
    
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(mfmatrix.getMatrix(pIndex,pIndex));

    PPDiagArray& ppDiag = ppMatrix.getDiag();
    PPDiagArray& ppCoeff = ppMatrix.getOffDiag();

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);

    const CRConnectivity& cr = ppMatrix.getConnectivity();

    const Array<int>& row = cr.getRow();

    ppDiag[0] = -1.;
    rCell[0]=0.;
    for(int nb=row[0]; nb<row[1]; nb++)
      ppCoeff[nb] = 0;
    
    if (mfmatrix.hasMatrix(pIndex,vIndex))
    {
        PVMatrix& pvMatrix =
          dynamic_cast<PVMatrix&>(mfmatrix.getMatrix(pIndex,vIndex));
        
        PVDiagArray& pvDiag = pvMatrix.getDiag();
        PVDiagArray& pvCoeff = pvMatrix.getOffDiag();

        pvDiag[0] = 0;
        for(int nb=row[0]; nb<row[1]; nb++)
          pvCoeff[nb] = 0;
    }
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
    }

    //    cout << "net boundary flux = " << netFlux << endl;
    
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
            const Array<int>& ibType = mesh.getIBType();
            const StorageSite& cells = mesh.getCells();
            const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
            for(int c=0; c<cells.getSelfCount(); c++)
              if (ibType[c] == Mesh::IBTYPE_FLUID)
                volumeSum += cellVolume[c];
        }

        netFlux /= volumeSum;

        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            const Array<int>& ibType = mesh.getIBType();
            const StorageSite& cells = mesh.getCells();
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
        
        _referencePP = pp[0];
    }
    else
      _referencePP = 0;
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
        
        if (_options.printNormalizedResiduals)
          cout << _niters << ": " << *mNormRatio << ";" << *cNormRatio <<  endl;
        else
          cout << _niters << ": " << *mNorm << ";" << *cNorm <<  endl;

        _niters++;
        if ((*mNormRatio < _options.momentumTolerance) &&
            (*cNormRatio < _options.continuityTolerance))
          return true;
    }
    return false;
  }

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
  
#ifndef USING_ATYPE_TANGENT
  
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
    const Array<int>& ibFaceIndices = mesh.getIBFaceList();
    const int nibf = ibFaces.getCount();
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
    const Array<int>& ibType = mesh.getIBType();
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
    
    if (nibf == 0)
      throw CException("getPressureIntegralonIBFaces: no IBFaces found!");
    return r;
  }

  VectorT3 getPVIntegral(const Field& velCoeffField,
                         const Mesh& mesh, const int faceGroupId)
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
            const VectorT3Array& velCoeff =
              dynamic_cast<const VectorT3Array&>(velCoeffField[faces]);
            const int nFaces = faces.getCount();
            const TArray& facePressure = dynamic_cast<const TArray&>(_flowFields.pressure[faces]);
            for(int f=0; f<nFaces; f++)
              r += velCoeff[f]*(faceArea[f]*facePressure[f]);

            found=true;
        }
    }
  if (!found)
    throw CException("getPVIntegral: invalid faceGroupID");
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

    const TArray& mu = dynamic_cast<const TArray&>(_flowFields.viscosity[cells]);

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
        
        stressTensor[n][0] = vgPlusTranspose[0][0]*mu[c] + pCell[c];
        stressTensor[n][1] = vgPlusTranspose[1][1]*mu[c] + pCell[c];
        stressTensor[n][2] = vgPlusTranspose[2][2]*mu[c] + pCell[c];
        stressTensor[n][3] = vgPlusTranspose[0][1]*mu[c];
        stressTensor[n][4] = vgPlusTranspose[1][2]*mu[c];
        stressTensor[n][5] = vgPlusTranspose[2][0]*mu[c];
    }

    return stressTensorPtr;
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
  T _referencePP;
  //AMG _momSolver;
  //AMG _continuitySolver;
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

template<class T>
bool
FlowModel<T>::advanceCoupled(const int niter)
{
  return _impl->advanceCoupled(niter);
}

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
FlowModel<T>::getPVIntegral(const Field& velCoeff,const Mesh& mesh, const int faceGroupId)
{
  return  _impl->getPVIntegral(velCoeff,mesh,faceGroupId);
}

template<class T>
Vector<T,3>
FlowModel<T>::getMomentumFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
 return  _impl->getMomentumFluxIntegral(mesh,faceGroupId);
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
FlowModel<T>::computeIBandSolidVelocity(const StorageSite& particles)
{
  return _impl->computeIBandSolidVelocity(particles);
}



template<class T>
boost::shared_ptr<ArrayBase>
FlowModel<T>::getStressTensor(const Mesh& mesh, const ArrayBase& cellIds)
{
  return _impl->getStressTensor(mesh, cellIds);
}


#ifndef USING_ATYPE_TANGENT
template<class T>  
void
FlowModel<T>::dumpContinuityMatrix(const string fileBase)
{
  _impl->dumpContinuityMatrix(fileBase);
}

#endif
