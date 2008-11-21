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

#include "Underrelaxer.h"
#include "MomentumPressureGradientDiscretization.h"
#include "AMG.h"
#include "Linearizer.h"
#include "GradientModel.h"

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


  typedef Array<Gradient<T> > PGradArray;

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

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            const FlowBC<T>& bc = *_bcMap[fg.id];

            VectorT3 bVelocity;
            bVelocity[0] = bc["specifiedXVelocity"];
            bVelocity[1] = bc["specifiedYVelocity"];
            bVelocity[2] = bc["specifiedZVelocity"];

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
                    massFlux[f] = density[c0]*dot(bVelocity,faceArea[f]);
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
    
            VectorT3 bVelocity;
            bVelocity[0] = bc["specifiedXVelocity"];
            bVelocity[1] = bc["specifiedYVelocity"];
            bVelocity[2] = bc["specifiedZVelocity"];

            if (bc.bcType == "NoSlipWall")
            {
                gbc.applyDirichletBC(bVelocity);
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
                        gbc.applyDirichletBC(f,bVelocity);
                    }
                }
            }
            else if ((bc.bcType == "Symmetry"))
            {
                VectorT3 zeroFlux(NumTypeTraits<VectorT3>::getZero());
                gbc.applyNeumannBC(zeroFlux);
            }
            else
              throw CException(bc.bcType + " not implemented for FlowModel");
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
    MFRPtr rNorm = _options.momentumLinearSolver->solve(ls);

    if (!_initialMomentumNorm) _initialMomentumNorm = rNorm;
        
    _options.momentumLinearSolver->cleanup();
    
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
    
    return rNorm;
  }


  void discretizeMassFluxInterior(const Mesh& mesh,
                                  const StorageSite& faces,
                                  MultiFieldMatrix& mfmatrix,
                                  const MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const TArray& faceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(mfmatrix.getMatrix(pIndex,pIndex));

    const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
    
    const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
    const VectorT3Array& Vprev = dynamic_cast<const VectorT3Array&>((*_previousVelocity)[cells]);

    const TArray& p = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
    const PGradArray& pGrad = dynamic_cast<const PGradArray&>(_flowFields.pressureGradient[cells]);

    const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);

    PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);
    PPDiagArray& ppDiag = ppMatrix.getDiag();

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

    const T momURF(_options["momentumURF"]);
    const T OneMinusmomURF(T(1.0)-momURF);
    
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
        const VectorT3& Af = faceArea[f];
        
        const T diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(Af,ds);

        const T momApBar0 = (momAp[c0][0]+momAp[c0][1]+momAp[c0][2])/3.0;
        const T momApBar1 = (momAp[c1][0]+momAp[c1][1]+momAp[c1][2])/3.0;
        const T momApBarFace = momApBar0 + momApBar1;

        const T VdotA0 = dot(V[c0],Af) - OneMinusmomURF*dot(Vprev[c0],Af);
        const T VdotA1 = dot(V[c1],Af) - OneMinusmomURF*dot(Vprev[c1],Af);

        
        const T dpf = cellVolume[c0]*(pGrad[c0]*ds) + cellVolume[c1]*(pGrad[c1]*ds);
        // const T dpf = (cellVolume[c0]*(pGrad[c0]*ds) + cellVolume[c1]*(pGrad[c1]*ds))/
        //(cellVolume[c0]+cellVolume[c1]) - p[c1] + p[c0];

        const T Vn = (VdotA0*momApBar0 + VdotA1*momApBar1 -dpf*diffMetric) / momApBarFace;
        //const T Vn = (VdotA0*momApBar0 + VdotA1*momApBar1 -
        //      (cellVolume[c0]+cellVolume[c1])*dpf*diffMetric) / momApBarFace;

        const T rhoF = 0.5*(rho[c0]+rho[c1]);
        const T aByMomAp = Af[0]*Af[0] / (momAp[c0][0] + momAp[c1][0]) +
          Af[1]*Af[1] / (momAp[c0][1] + momAp[c1][1]) +
          Af[2]*Af[2] / (momAp[c0][2] + momAp[c1][2]);


        const T pCoeff = rhoF*aByMomAp*(cellVolume[c0]+cellVolume[c1])/(dot(Af,ds));

        massFlux[f] = rhoF*Vn - pCoeff*(p[c0]-p[c1]) + (1-momURF)*massFlux[f];
        //massFlux[f] = rhoF*Vn  + OneMinusmomURF*massFlux[f];

        rCell[c0] -= massFlux[f];
        rCell[c1] += massFlux[f];
        
        ppAssembler.getCoeff01(f) -=pCoeff;
        ppAssembler.getCoeff10(f) -=pCoeff;

        ppDiag[c0] += pCoeff;
        ppDiag[c1] += pCoeff;
    }

#if 1
    if (mfmatrix.hasMatrix(pIndex,vIndex))
    {
        PVMatrix& pvMatrix =
          dynamic_cast<PVMatrix&>(mfmatrix.getMatrix(pIndex,vIndex));
        
        PVAssembler& pvAssembler = pvMatrix.getPairWiseAssembler(faceCells);
        PVDiagArray& pvDiag = pvMatrix.getDiag();

        for(int f=0; f<nFaces; f++)
        {
            const int c0 = faceCells(f,0);
            const int c1 = faceCells(f,1);
            const VectorT3& Af = faceArea[f];
            
            const T momApBar0 = (momAp[c0][0]+momAp[c0][1]+momAp[c0][2])/3.0;
            const T momApBar1 = (momAp[c1][0]+momAp[c1][1]+momAp[c1][2])/3.0;
            const T momApBarFace = momApBar0 + momApBar1;
            const T rhoF = 0.5*(rho[c0]+rho[c1]);

            VectorT3T coeff0(rhoF*momApBar0/momApBarFace*Af);
            VectorT3T coeff1(rhoF*momApBar1/momApBarFace*Af);
            
        
            pvAssembler.getCoeff01(f) -=coeff1;
            pvAssembler.getCoeff10(f) +=coeff0;
            
            pvDiag[c0] -= coeff0;
            pvDiag[c1] += coeff1;
        }
    }
#endif
  }

  T fixedFluxContinuityBC(const StorageSite& faces,
                          const Mesh& mesh,
                          MultiFieldMatrix& matrix,
                          MultiField& xField,
                          MultiField& rField,
                          const FlowBC<T>& bc)
  {
    const StorageSite& cells = mesh.getCells();

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(matrix.getMatrix(pIndex,pIndex));

    FMatrix& dFluxdP = dynamic_cast<FMatrix&>(matrix.getMatrix(mfIndex,pIndex));

    PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);
    PPDiagArray& ppDiag = ppMatrix.getDiag();

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
    const TArray& density = dynamic_cast<const TArray&>(_flowFields.density[cells]);

    VectorT3 bVelocity;
    bVelocity[0] = bc["specifiedXVelocity"];
    bVelocity[1] = bc["specifiedYVelocity"];
    bVelocity[2] = bc["specifiedZVelocity"];

    const int nFaces = faces.getCount();

    T netFlux(0.);
    
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        massFlux[f] = density[c0]*dot(bVelocity,faceArea[f]);

        rCell[c0] -= massFlux[f];

        netFlux += massFlux[f];
        ppAssembler.getCoeff01(f) =0;
        ppAssembler.getCoeff10(f) =1;
        ppDiag[c1] = -1;
        rCell[c1] = 0.;
        ppMatrix.setBoundary(c1);

        dFluxdP.setCoeffL(f,T(0.));
        dFluxdP.setCoeffR(f,T(0.));
        
    }
#if 1
    if (matrix.hasMatrix(vIndex,pIndex))
    {
        VPMatrix& vpMatrix =
          dynamic_cast<VPMatrix&>(matrix.getMatrix(vIndex,pIndex));
        
        VPAssembler& vpAssembler = vpMatrix.getPairWiseAssembler(faceCells);
        VPDiagArray& vpDiag = vpMatrix.getDiag();

        for(int f=0; f<nFaces; f++)
        {
            const int c0 = faceCells(f,0);
            vpDiag[c0] += vpAssembler.getCoeff01(f);
            vpAssembler.getCoeff01(f) = 0;
        }
    }
#endif
    return netFlux;
  }

  T fixedPressureContinuityBC(const StorageSite& faces,
                              const Mesh& mesh,
                              MultiFieldMatrix& matrix,
                              const MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    MultiField::ArrayIndex mfIndex(&_flowFields.massFlux,&faces);

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(matrix.getMatrix(pIndex,pIndex));

    const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
    
    const VectorT3Array& V = dynamic_cast<const VectorT3Array&>(_flowFields.velocity[cells]);
    const VectorT3Array& Vprev = dynamic_cast<const VectorT3Array&>((*_previousVelocity)[cells]);

    const TArray& p = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
    const PGradArray& pGrad = dynamic_cast<const PGradArray&>(_flowFields.pressureGradient[cells]);

    const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);

    PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);
    PPDiagArray& ppDiag = ppMatrix.getDiag();

    TArray& rCell = dynamic_cast<TArray&>(rField[pIndex]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

    FMatrix& dFluxdP = dynamic_cast<FMatrix&>(matrix.getMatrix(mfIndex,pIndex));

    const T momURF(_options["momentumURF"]);
    const T OneMinusmomURF(T(1.0)-momURF);
    
    const int nFaces = faces.getCount();

    T netFlux(0.);
    
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
        const VectorT3& Af = faceArea[f];
        
        const T dpf = pGrad[c0]*ds - p[c1] + p[c0] ;
        const T rhoF = rho[c0];

        // Q < 0
        const T Q = rhoF*(Af[0]*Af[0] / momAp[c0][0]  +
                                 Af[1]*Af[1] / momAp[c0][1]  +
                                 Af[2]*Af[2] / momAp[c0][2])*cellVolume[c0]/(dot(Af,ds));
        
        const T massFluxI = rhoF*(dot(V[c0],Af) - OneMinusmomURF*dot(Vprev[c0],Af)) - Q*dpf +
          OneMinusmomURF*massFlux[f];

        const VectorT3& Vb = V[c1];
        const T massFluxB = rhoF*dot(Vb,Af);
        

        T Vb_dpdVb(0);
        if (massFluxB < 0)
          Vb_dpdVb = mag2(Vb)*rhoF;

        netFlux += massFluxI;
        
        massFlux[f] = massFluxI;

        dFluxdP.setCoeffL(f,Q);
        dFluxdP.setCoeffR(f,-Q);

        // contribution to cell equation
        rCell[c0] -= massFlux[f];
        ppDiag[c0] += Q;
        ppAssembler.getCoeff01(f) = 0.;

        if (massFluxB - Vb_dpdVb*Q != 0.0)
        {
            // equation for boundary 
            ppAssembler.getCoeff10(f) = -Vb_dpdVb * Q / (massFluxB - Vb_dpdVb*Q);
            rCell[c1] = Vb_dpdVb*(massFluxI-massFluxB) / (massFluxB - Vb_dpdVb*Q);
            rCell[c1] =0;
            ppDiag[c1] = -1;
            
            // eliminate boundary dependence from cell equation
            ppAssembler.getCoeff01(f) = 0;
            if (ppAssembler.getCoeff10(f) > 0)
              ppDiag[c0] += ppAssembler.getCoeff10(f)*Q;
            rCell[c0] += rCell[c1]*Q;

        }
        else
        {
            // treat as fixed pressure
            dFluxdP.setCoeffR(f,T(0.));
            ppDiag[c1] = -1;
            rCell[c1] = 0;
            ppAssembler.getCoeff10(f) = 0.;
            ppAssembler.getCoeff01(f) = 0.;
        }
        ppMatrix.setBoundary(c1);
    }
    return netFlux;
  }

  void pressureBoundaryPostContinuitySolve(const StorageSite& faces,
                                           const Mesh& mesh,
                                           const T bp)
  {
    const StorageSite& cells = mesh.getCells();

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    
    const TArray& faceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
    

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
    VectorT3Array& V = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
    TArray& p = dynamic_cast<TArray&>(_flowFields.pressure[cells]);
    const TArray& rho = dynamic_cast<TArray&>(_flowFields.density[cells]);
    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        const T rhoF = rho[c0];

        const T Vn = massFlux[f]/(rhoF*faceAreaMag[f]);

        if (massFlux[f] > 0)
        {
            p[c1]=bp;
        }
        else
        {
            V[c1] = -Vn*faceArea[f]/faceAreaMag[f];
            p[c1] = bp - 0.5*rhoF*mag2(V[c1]);
        }
    }
  }

  void correctVelocityInterior(const Mesh& mesh,
                               const StorageSite& faces,
                               const MultiField& ppField)                               
  {
    const StorageSite& cells = mesh.getCells();

    const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    const VectorT3Array& cellCentroid =  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
    
    VectorT3Array& V = dynamic_cast<VectorT3Array&>(_flowFields.velocity[cells]);
    const TArray& pp = dynamic_cast<const TArray&>(ppField[pIndex]);
    const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
        const VectorT3& Af = faceArea[f];

        const T aByMomAp0 = Af[0]*Af[0] / momAp[c0][0] +
          Af[1]*Af[1] / momAp[c0][1] +
          Af[2]*Af[2] / momAp[c0][2];

        const T aByMomAp1 = Af[0]*Af[0] / momAp[c1][0] +
          Af[1]*Af[1] / momAp[c1][1] +
          Af[2]*Af[2] / momAp[c1][2];

        const T Adotes = dot(Af,ds)/mag(ds);
        const T coeff0  = cellVolume[c0]*rho[c0]*aByMomAp0/Adotes;
        const T coeff1  = cellVolume[c1]*rho[c1]*aByMomAp1/Adotes;
        
        const T ppFace = (coeff0*pp[c0]+coeff1*pp[c1])/(coeff0+coeff1);
        const VectorT3 ppA = ppFace*faceArea[f];

        V[c0] += ppA/momAp[c0];
        V[c1] -= ppA/momAp[c1];
    }
  }

  void updateFacePressureInterior(const Mesh& mesh,
                                  const StorageSite& faces)
  {
    const StorageSite& cells = mesh.getCells();

    const VectorT3Array& faceArea = dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
    const VectorT3Array& cellCentroid =  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    const VVDiagArray& momAp = dynamic_cast<const VVDiagArray&>((*_momApField)[cells]);
    
    const TArray& pCell = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
    TArray& pFace = dynamic_cast<TArray&>(_flowFields.pressure[faces]);
    const TArray& rho = dynamic_cast<const TArray&>(_flowFields.density[cells]);
    const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        const VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
        const VectorT3& Af = faceArea[f];

        const T aByMomAp0 = Af[0]*Af[0] / momAp[c0][0] +
          Af[1]*Af[1] / momAp[c0][1] +
          Af[2]*Af[2] / momAp[c0][2];

        const T aByMomAp1 = Af[0]*Af[0] / momAp[c1][0] +
          Af[1]*Af[1] / momAp[c1][1] +
          Af[2]*Af[2] / momAp[c1][2];

        const T Adotes = dot(Af,ds)/mag(ds);
        const T coeff0  = cellVolume[c0]*rho[c0]*aByMomAp0/Adotes;
        const T coeff1  = cellVolume[c1]*rho[c1]*aByMomAp1/Adotes;
        
        pFace[f] = (coeff0*pCell[c0]+coeff1*pCell[c1])/(coeff0+coeff1);
    }
  }

  // the real face pressure is calculated by bc's and stored in p[c1]
  // for boundary faces but if we are also storing the pressure at all
  // faces we need to copy it here
  
  void updateFacePressureBoundary(const Mesh& mesh,
                                  const StorageSite& faces)
  {
    const StorageSite& cells = mesh.getCells();

    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    const TArray& pCell = dynamic_cast<const TArray&>(_flowFields.pressure[cells]);
    TArray& pFace = dynamic_cast<TArray&>(_flowFields.pressure[faces]);

    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        //        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        
        pFace[f] = pCell[c1];
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

  void correctMassFluxInterior(const Mesh& mesh,
                               const StorageSite& faces,
                               MultiFieldMatrix& mfmatrix,
                               const MultiField& xField)
  {
    const StorageSite& cells = mesh.getCells();
    const CRConnectivity& faceCells = mesh.getFaceCells(faces);

    MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
    const TArray& pp = dynamic_cast<const TArray&>(xField[pIndex]);

    PPMatrix& ppMatrix =
      dynamic_cast<PPMatrix&>(mfmatrix.getMatrix(pIndex,pIndex));
    
    PPAssembler& ppAssembler = ppMatrix.getPairWiseAssembler(faceCells);

    TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
    
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        massFlux[f] -= ppAssembler.getCoeff01(f)*pp[c1] -
          ppAssembler.getCoeff10(f)*pp[c0];
    }

#if 1
    MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    if (mfmatrix.hasMatrix(pIndex,vIndex))
    {
        PVMatrix& pvMatrix =
          dynamic_cast<PVMatrix&>(mfmatrix.getMatrix(pIndex,vIndex));
        
        PVAssembler& pvAssembler = pvMatrix.getPairWiseAssembler(faceCells);
        const VectorT3Array& Vp = dynamic_cast<const VectorT3Array&>(xField[vIndex]);

        for(int f=0; f<nFaces; f++)
        {
            const int c0 = faceCells(f,0);
            const int c1 = faceCells(f,1);
            
            massFlux[f] += pvAssembler.getCoeff01(f)*Vp[c1] +
              pvAssembler.getCoeff10(f)*Vp[c0];
        }
    }
#endif
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
    
#if 0
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        _flowFields.pressureGradient.syncGather(cells);
    }
#endif
    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const StorageSite& iFaces = mesh.getInteriorFaceGroup().site;
        
        MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
        
        // interior
        
        discretizeMassFluxInterior(mesh,iFaces,matrix,x,b);
        
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

    if (this->_useReferencePressure)
    {
        T volumeSum(0.);
        
        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            const StorageSite& cells = mesh.getCells();
            const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
            for(int c=0; c<cells.getSelfCount(); c++) volumeSum += cellVolume[c];
        }

        netFlux /= volumeSum;

        for (int n=0; n<numMeshes; n++)
        {
            const Mesh& mesh = *_meshes[n];
            const StorageSite& cells = mesh.getCells();
            const TArray& cellVolume = dynamic_cast<const TArray&>(_geomFields.volume[cells]);
            
            MultiField::ArrayIndex pIndex(&_flowFields.pressure,&cells);
            TArray& rCell = dynamic_cast<TArray&>(b[pIndex]);
            
            for(int c=0; c<cells.getSelfCount(); c++)
            {
                rCell[c] += netFlux*cellVolume[c];
            }
        }
    }
    setDirichlet(matrix,b);
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
        else
          correctVelocityInterior(mesh,iFaces,ppField);
        
        updateFacePressureInterior(mesh,iFaces);
        
        // boundary
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            const FlowBC<T>& bc = *_bcMap[fg.id];
              
            correctMassFluxBoundary(faces,ppField);
            
            if (!coupled)
              correctVelocityBoundary(mesh,faces,ppField);

            if (bc.bcType == "PressureBoundary")
            {
                T bp = bc["specifiedPressure"];
                pressureBoundaryPostContinuitySolve(faces,mesh,bp);
            }
            
            updateFacePressureBoundary(mesh,faces);
        }

    }

#if 0
    for(int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        
        const StorageSite& cells = mesh.getCells();
        //        MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
        // const bool coupled =  matrix.hasMatrix(pIndex,vIndex);
        //if (!coupled)
        _flowFields.velocity.syncGather(cells);
    }
#endif
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

    MFRPtr rNorm = _options.pressureLinearSolver->solve(*ls);

    if (!_initialContinuityNorm) _initialContinuityNorm = rNorm;
        
    ls->postSolve();
    _options.pressureLinearSolver->cleanup();
    
    postContinuitySolve(*ls);

    // discard the momentum ap coeffficients
    
    _momApField = shared_ptr<Field>();
    return rNorm;
  }


  void advance(const int niter)
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
          break;
    }
  }

  void advanceCoupled(const int niter)
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
          break;
    }
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
            cout << "   " << vp.first << " "  << vp.second <<  endl;
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
  
  //  LinearSolver& getMomentumSolver() {return _momSolver;}
  // LinearSolver& getContinuitySolver() {return _continuitySolver;}
  
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
void
FlowModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}

template<class T>
void
FlowModel<T>::advanceCoupled(const int niter)
{
  _impl->advanceCoupled(niter);
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
FlowModel<T>::getMomentumFluxIntegral(const Mesh& mesh, const int faceGroupId)
{
 return  _impl->getMomentumFluxIntegral(mesh,faceGroupId);
}

#ifndef USING_ATYPE_TANGENT
template<class T>  
void
FlowModel<T>::dumpContinuityMatrix(const string fileBase)
{
  _impl->dumpContinuityMatrix(fileBase);
}

#endif
