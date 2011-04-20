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
//#include "SourceDiscretization.h"
#include "SourceDiscretizationene.h"
#include "SourceDiscretizationdissi.h"

template<class T>
class  KeModel<T>::Impl
{
public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Gradient<VectorT3> VGradType;
  typedef Array<Gradient<VectorT3> > VGradArray;
  typedef Gradient<T> EGradType;
  typedef Array<EGradType> EGradArray;
  typedef Gradient<T> DGradType;
  typedef Array<DGradType> DGradArray;
 
  typedef CRMatrix<T,T,T> T_Matrix;
  
  Impl(const GeomFields& geomFields,
       KeFields& keFields,
       FlowFields& flowFields,
       const MeshList& meshes) :
    _meshes(meshes),
    _geomFields(geomFields),
    _keFields(keFields),
    _flowFields(flowFields),
    _velocityGradientModel(_meshes,_flowFields.velocity,
                           _flowFields.velocityGradient,_geomFields),
    _energyGradientModel(_meshes,_keFields.energy,
                         _keFields.energyGradient,_geomFields),
    _dissipationGradientModel(_meshes,_keFields.dissipation,
                              _keFields.dissipationGradient,_geomFields),
    _initialNormk(),
    _initialNorm(),
    _niters(0)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
  /*      
	FlowVC<T> *vc1(new FlowVC<T>());
        vc1->vcType = "flow";
        _vcMap[mesh.getID()] = vc1;
*/
        KeVC<T> *vc(new KeVC<T>());
        vc->vcType = "flow";
	_vcMap[mesh.getID()] = vc;
        
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            KeBC<T> *bc(new KeBC<T>());
            
            _bcMap[fg.id] = bc;
	    
            if ((fg.groupType == "symmetry"))
	      {
		bc->bcType = "Symmetry";
	      }

            else if ((fg.groupType == "velocity-inlet"))
            {
              bc->bcType = "Specifiedkandepsilon";
            }
             

             else if ((fg.groupType == "pressure-inlet") ||
                      (fg.groupType == "pressure-outlet"))
              {
                bc->bcType = "PressureBoundary";
              }

            else
	      throw CException("KeModel: unknown face group type "
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
        //const FlowVC<T>&vc1 = *_vcMap[mesh.getID()];
        const KeVC<T>&vc = *_vcMap[mesh.getID()];

        const StorageSite& cells = mesh.getCells();
        //energy(k)

        shared_ptr<TArray> kCell(new TArray(cells.getCount()));
        *kCell = vc["InitialEnergy"];
        _keFields.energy.addArray(cells,kCell);

	//dissipation(e)
        shared_ptr<TArray> eCell(new TArray(cells.getCount()));
        *eCell = vc["InitialDissipation"];
        _keFields.dissipation.addArray(cells,eCell);

    //    if (_options.transient)
  //      {
            _keFields.dissipationN1.addArray(cells,
                                            dynamic_pointer_cast<ArrayBase>(eCell->newCopy()));
            if (_options.timeDiscretizationOrder > 1)
              _keFields.dissipationN2.addArray(cells,
                                              dynamic_pointer_cast<ArrayBase>(eCell->newCopy()));

//        }


	
	//initial velocity gradient array
	shared_ptr<VGradArray> gradV(new VGradArray(cells.getCount()));
	gradV->zero();
	_flowFields.velocityGradient.addArray(cells,gradV);
        
	 //initial energy gradient array
        shared_ptr<EGradArray> gradE(new EGradArray(cells.getCount()));
        gradE->zero();
        _keFields.energyGradient.addArray(cells,gradE);

	 //initial dissipation gradient array
        shared_ptr<DGradArray> gradD(new DGradArray(cells.getCount()));
        gradD->zero();
        _keFields.dissipationGradient.addArray(cells,gradD);
       
       // density field
        shared_ptr<TArray> densityCell(new TArray(cells.getCount()));
        *densityCell = T(1.225);
        //*densityCell = vc1["density"];
        _flowFields.density.addArray(cells,densityCell);
 
        //c1
        shared_ptr<TArray> c1Cell(new TArray(cells.getCount()));
        *c1Cell = vc["c1"];
        _keFields.c1.addArray(cells,c1Cell);

        //c2
        shared_ptr<TArray> c2Cell(new TArray(cells.getCount()));
        *c2Cell = vc["c2"];
        _keFields.c2.addArray(cells,c2Cell);
        
        //sigmak
        shared_ptr<TArray> sigmakCell(new TArray(cells.getCount()));
        *sigmakCell = vc["sigmak"];
        _keFields.sigmak.addArray(cells,sigmakCell);

        //sigmae
        shared_ptr<TArray> sigmaeCell(new TArray(cells.getCount()));
        *sigmaeCell = vc["sigmae"];
        _keFields.sigmae.addArray(cells,sigmaeCell);

        //cmu
        shared_ptr<TArray> cmuCell(new TArray(cells.getCount()));
        *cmuCell = vc["cmu"];
        _keFields.cmu.addArray(cells,cmuCell);


        //kflux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _keFields.kFlux.addArray(faces,fluxFace);

        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _keFields.kFlux.addArray(faces,fluxFace);

        }
    
        //eflux at faces
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _keFields.eFlux.addArray(faces,fluxFace);

        }
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            shared_ptr<TArray> fluxFace(new TArray(faces.getCount()));

            fluxFace->zero();
            _keFields.eFlux.addArray(faces,fluxFace);

        }
 

       //laminar viscosity
        shared_ptr<TArray> lmuCell(new TArray(cells.getCount()));
        *lmuCell = T(1.7894e-5);
       // *lmuCell = vc["viscosity"];
        _flowFields.viscosity.addArray(cells,lmuCell);

       //turbulent viscosity
        shared_ptr<TArray> muCell(new TArray(cells.getCount()));
       *muCell = T(1e-05); 
       //*muCell = vc["eddyviscosity"];
        _flowFields.eddyviscosity.addArray(cells,muCell);


        //total viscosity
        shared_ptr<TArray> tmuCell(new TArray(cells.getCount()));
        *tmuCell = T(1e-03);
        // *tmuCell = vc["totalviscosity"];
        _flowFields.totalviscosity.addArray(cells,tmuCell);

    }
    //_keFields.energy.synclocal();
   // _keFields.dissipation.synclocal(); 
    _niters  =0;
    _initialNormk = MFRPtr();
    _initialNorm = MFRPtr();

  }
  
  KeBCMap& getBCMap() {return _bcMap;}
  KeVCMap& getVCMap() {return _vcMap;}
 // FlowVCMap& getVCMap() {return _vcMap;}
  KeBC<T>& getBC(const int id) {return *_bcMap[id];}

  KeModelOptions<T>& getOptions() {return _options;}

  void updateTimek()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        TArray& k =
          dynamic_cast<TArray&>(_keFields.energy[cells]);
        TArray& kN1 =
          dynamic_cast<TArray&>(_keFields.energyN1[cells]);

        if (_options.timeDiscretizationOrder > 1)
        {
            TArray& kN2 =
              dynamic_cast<TArray&>(_keFields.energyN2[cells]);
            kN2 = kN1;
        }
        kN1 = k;

    }
  }




  void updateTimee()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        TArray& e =
          dynamic_cast<TArray&>(_keFields.dissipation[cells]);
        TArray& eN1 =
          dynamic_cast<TArray&>(_keFields.dissipationN1[cells]);

        if (_options.timeDiscretizationOrder > 1)
        {
            TArray& eN2 =
              dynamic_cast<TArray&>(_keFields.dissipationN2[cells]);
            eN2 = eN1;
        }
        eN1 = e;
    }
  }




  void initLinearizationk(LinearSystem& lsk)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex kIndex(&_keFields.energy,&cells);

        lsk.getX().addArray(kIndex,_keFields.energy.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        lsk.getMatrix().addMatrix(kIndex,kIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_keFields.kFlux,&faces);
            lsk.getX().addArray(fIndex,_keFields.kFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            lsk.getMatrix().addMatrix(fIndex,kIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            lsk.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_keFields.kFlux,&faces);
            lsk.getX().addArray(fIndex,_keFields.kFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            lsk.getMatrix().addMatrix(fIndex,kIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            lsk.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

    }
  }

  void linearizeenergy(LinearSystem& lsk)
  {
    _velocityGradientModel.compute();
    _energyGradientModel.compute();
    DiscrList discretizations;

    if (_options.transient)
    {

       shared_ptr<Discretization>
         td(new TimeDerivativeDiscretization<T,T,T>
           (_meshes,_geomFields,
            _keFields.energy,
            _keFields.energyN1,
            _keFields.energyN2,
            _flowFields.density,
            _options["timeStep"]));

       discretizations.push_back(td);
    }   
    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
	 (_meshes,_geomFields,
	  _keFields.energy,
          _keFields.c1,
	  _keFields.energyGradient));
    discretizations.push_back(dd);
    
    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<T,T,T>
	 (_meshes,_geomFields,
          _keFields.energy,
	  _flowFields.massFlux,
          _flowFields.continuityResidual,
	  _keFields.energyGradient,
          _options.useCentralDifference));
    discretizations.push_back(cd);
    
    
    shared_ptr<Discretization>
      sd(new SourceDiscretizationene<T>
	 (_meshes, 
	  _geomFields,
          _keFields.energy, 
          _flowFields.velocity,
	  _flowFields.eddyviscosity,
          _keFields.dissipation,
          _flowFields.velocityGradient));
    discretizations.push_back(sd);
   
    shared_ptr<Discretization>
      ene(new GenericIBDiscretization<T,T,T>
	  (_meshes,_geomFields,_keFields.energy));
    
    discretizations.push_back(ene);
    
    Linearizer linearizer;
    
    linearizer.linearize(discretizations,_meshes,lsk.getMatrix(),
                         lsk.getX(), lsk.getB());
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
	
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
           // const int nFaces = faces.getCount();

            const KeBC<T>& bc = *_bcMap[fg.id];
            

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _keFields.energy,
                                  _keFields.kFlux,
                                  lsk.getMatrix(), lsk.getX(), lsk.getB());

            FloatValEvaluator<T>
             bk(bc.getVal("Specifiedk"),faces);
            TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);

            const int nFaces = faces.getCount();

	    
          //  if (bc.bcType == "Specifiedk")
          if ( bc.bcType == "Specifiedkandepsilon")
	      {
                for(int f=0; f<nFaces; f++)
		   {
                     if (massFlux[f] > 0.)
	               {
                         gbc.applyExtrapolationBC(f);
		       }
                      else
			{
                          gbc.applyDirichletBC(f,bk[f]);
			}
		     }
		
               }
           else if ((bc.bcType == "PressureBoundary"))
            {
            for(int f=0; f<nFaces; f++)
                {
                    if (massFlux[f] > 0.)
                    {
                        gbc.applyExtrapolationBC(f);
                    }
                    else
                    {
                        gbc.applyDirichletBC(f,bk[f]);
                    }
                }
              }
 



           else if ((bc.bcType == "Symmetry"))
	      {
                 gbc.applySymmetryBC();
	      }
            else
              throw CException(bc.bcType + " not implemented for KeModel");
	  }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _keFields.energy,
                                  _keFields.kFlux,
                                  lsk.getMatrix(), lsk.getX(), lsk.getB());

            gbc.applyInterfaceBC();

	}
      }
      DiscrList discretizations2;
      shared_ptr<Discretization>
        ud(new Underrelaxer<T,T,T>
         (_meshes,_keFields.energy,
          _options["energyURF"]));

    discretizations2.push_back(ud);

    linearizer.linearize(discretizations2,_meshes,lsk.getMatrix(),
                         lsk.getX(), lsk.getB());
 
  }
 
  void initLinearization(LinearSystem& lse)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();
        MultiField::ArrayIndex eIndex(&_keFields.dissipation,&cells);

        lse.getX().addArray(eIndex,_keFields.dissipation.getArrayPtr(cells));

        const CRConnectivity& cellCells = mesh.getCellCells();

        shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));

        lse.getMatrix().addMatrix(eIndex,eIndex,m);

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_keFields.eFlux,&faces);
            lse.getX().addArray(fIndex,_keFields.eFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            lse.getMatrix().addMatrix(fIndex,eIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            lse.getMatrix().addMatrix(fIndex,fIndex,mff);
        }
         foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;

            MultiField::ArrayIndex fIndex(&_keFields.eFlux,&faces);
            lse.getX().addArray(fIndex,_keFields.eFlux.getArrayPtr(faces));

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);

            shared_ptr<Matrix> mft(new FluxJacobianMatrix<T,T>(faceCells));
            lse.getMatrix().addMatrix(fIndex,eIndex,mft);

            shared_ptr<Matrix> mff(new DiagonalMatrix<T,T>(faces.getCount()));
            lse.getMatrix().addMatrix(fIndex,fIndex,mff);
        }

    }
  }

 
  void linearizedissipation(LinearSystem& lse)
  {
    _velocityGradientModel.compute();
    _dissipationGradientModel.compute();

    DiscrList discretizations1;

    shared_ptr<Discretization>
      td(new TimeDerivativeDiscretization<T,T,T>
        (_meshes,_geomFields,
         _keFields.dissipation,
         _keFields.dissipationN1,
         _keFields.dissipationN2,
         _flowFields.density,
         _options["timeStep"]));
    discretizations1.push_back(td);
    

    shared_ptr<Discretization>
      dd(new DiffusionDiscretization<T,T,T>
         (_meshes,_geomFields,
          _keFields.dissipation,
          _keFields.c2,
          _keFields.dissipationGradient));
    discretizations1.push_back(dd);

    shared_ptr<Discretization>
      cd(new ConvectionDiscretization<T,T,T>
         (_meshes,_geomFields,
          _keFields.dissipation,
          _flowFields.massFlux,
          _flowFields.continuityResidual,
          _keFields.dissipationGradient,
          _options.useCentralDifference));
    discretizations1.push_back(cd);


    shared_ptr<Discretization>
      sd(new SourceDiscretizationdissi<T,T,T>
         (_meshes,
          _geomFields,
          _keFields.dissipation,
          _flowFields.velocity,
          _flowFields.eddyviscosity,
          _keFields.energy,
          _flowFields.velocityGradient));
    discretizations1.push_back(sd);

    shared_ptr<Discretization>
      dis(new GenericIBDiscretization<T,T,T>
          (_meshes,_geomFields,_keFields.dissipation));

    discretizations1.push_back(dis);

    Linearizer linearizer;
    linearizer.linearize(discretizations1,_meshes,lse.getMatrix(),
                         lse.getX(), lse.getB());

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
           // const int nFaces = faces.getCount();

            const KeBC<T>& bc = *_bcMap[fg.id];


            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _keFields.dissipation,
                                  _keFields.eFlux,
                                  lse.getMatrix(), lse.getX(), lse.getB());

            FloatValEvaluator<T>
              be(bc.getVal("specifiede"),faces);
           
            TArray& massFlux = dynamic_cast<TArray&>(_flowFields.massFlux[faces]);
            const int nFaces = faces.getCount();

           // if (bc.bcType == "Specifiede")
            if (bc.bcType == "Specifiedkandepsilon")
            {
              for(int f=0; f<nFaces; f++)
               {
                if (massFlux[f] > 0.)
                 {
                   gbc.applyExtrapolationBC(f);
                 }
                else
                 {
                   gbc.applyDirichletBC(f,be[f]);
                 }
               }
                                
            }

            else if ((bc.bcType == "PressureBoundary"))
            {
            for(int f=0; f<nFaces; f++)
                {
                    if (massFlux[f] > 0.)
                    {
                        gbc.applyExtrapolationBC(f);
                    }
                    else
                    {
                        gbc.applyDirichletBC(f,be[f]);
                    }
                }
              }

            else if ((bc.bcType == "Symmetry"))
            {
                   gbc.applySymmetryBC();
            }
            else
              throw CException(bc.bcType + " not implemented for KeModel");
        }

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,
                                  _keFields.dissipation,
                                  _keFields.eFlux,
                                  lse.getMatrix(), lse.getX(), lse.getB());

            gbc.applyInterfaceBC();

 
       }
    }
      DiscrList discretizations3;
      shared_ptr<Discretization>
        ud1(new Underrelaxer<T,T,T>
         (_meshes,_keFields.dissipation,
          _options["dissipationURF"]));

    discretizations3.push_back(ud1);

    linearizer.linearize(discretizations3,_meshes,lse.getMatrix(),
                         lse.getX(), lse.getB());
 
}


  void getViscosity(const Mesh& mesh)
  {
    const StorageSite& cells = mesh.getCells();

//turbulent viscosity
     TArray& muCell =
      dynamic_cast<TArray&>(_flowFields.eddyviscosity[cells]);

//laminar viscosity
     TArray& lmuCell =
      dynamic_cast<TArray&>(_flowFields.viscosity[cells]);

//total viscosity
     TArray& tmuCell =
      dynamic_cast<TArray&>(_flowFields.totalviscosity[cells]);

    const TArray& eCell =
      dynamic_cast<const TArray&>(_keFields.dissipation[cells]);

    const TArray& kCell =
      dynamic_cast<const TArray&>(_keFields.energy[cells]);

    TArray& c1Cell =
     dynamic_cast<TArray&>(_keFields.c1[cells]);

    TArray& c2Cell =
      dynamic_cast<TArray&>(_keFields.c2[cells]);

    const TArray& sigmakCell =
      dynamic_cast<const TArray&>(_keFields.sigmak[cells]);

    const TArray& sigmaeCell =
      dynamic_cast<const TArray&>(_keFields.sigmae[cells]);

    const TArray& cmuCell =
      dynamic_cast<const TArray&>(_keFields.cmu[cells]);


    const int nCells = cells.getCount();

    for(int c=0; c<nCells; c++)
    {
        muCell[c] = (cmuCell[c]*kCell[c]*kCell[c])/eCell[c];
        c1Cell[c] = muCell[c]/sigmakCell[c];
        c2Cell[c] = muCell[c]/sigmaeCell[c];
        tmuCell[c] = muCell[c] + lmuCell[c];
    }

  }
   

  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)
    { 
      {
        LinearSystem lsk;
        initLinearizationk(lsk);
        
        lsk.initAssembly();

        linearizeenergy(lsk);

        lsk.initSolve();

        MFRPtr rNorm(_options.getLinearSolver().solve(lsk));

        if (!_initialNormk) _initialNormk = rNorm;
        
        MFRPtr normRatio((*rNorm)/(*_initialNormk));

        cout << _niters << ": " << *rNorm << endl;

        
        _options.getLinearSolver().cleanup();

        lsk.postSolve();
        lsk.updateSolution();

        _niters++;
        if (*rNorm < _options.absoluteTolerance ||
            *normRatio < _options.relativeTolerance)
          break;
      }
      { 
        LinearSystem lse;
        initLinearization(lse);

        lse.initAssembly();

        linearizedissipation(lse);

        lse.initSolve();

        MFRPtr rNorm(_options.getLinearSolver().solve(lse));

        if (!_initialNorm) _initialNorm = rNorm;

        MFRPtr normRatio((*rNorm)/(*_initialNorm));

        cout << _niters << ": " << *rNorm << endl;


        _options.getLinearSolver().cleanup();

        lse.postSolve();
        lse.updateSolution();

        _niters++;
        if (*rNorm < _options.absoluteTolerance ||
            *normRatio < _options.relativeTolerance)
          break;
    }
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        getViscosity(mesh);
    }
    }
  }
  
 void printBCs()
  {
    foreach(typename KeBCMap::value_type& pos, _bcMap)
    {
        cout << "Face Group " << pos.first << ":" << endl;
        cout << "    bc type " << pos.second->bcType << endl;
        foreach(typename KeBC<T>::value_type& vp, *pos.second)
        {
            cout << "   " << vp.first << " "  << vp.second.constant <<  endl;
        }
    }
  }

private:
  const MeshList _meshes;
  const GeomFields& _geomFields;
  KeFields& _keFields;
  FlowFields& _flowFields;

  KeBCMap _bcMap;
  
  KeVCMap _vcMap;
  KeModelOptions<T> _options;
  GradientModel<VectorT3> _velocityGradientModel;  
  GradientModel<T> _energyGradientModel;
  GradientModel<T> _dissipationGradientModel;
  MFRPtr _initialNormk; 
  MFRPtr _initialNorm;
  int _niters;
};

template<class T>
KeModel<T>::KeModel(const GeomFields& geomFields,
                    KeFields& keFields,
                    FlowFields& flowFields,
                    const MeshList& meshes) :
  Model(meshes),
  _impl(new Impl(geomFields,keFields,flowFields,meshes))
{
  logCtor();
}


template<class T>
KeModel<T>::~KeModel()
{
  logDtor();
}

template<class T>
void
KeModel<T>::init()
{
  _impl->init();
}
  
template<class T>
typename KeModel<T>::KeBCMap&
KeModel<T>::getBCMap() {return _impl->getBCMap();}

template<class T>
typename KeModel<T>::KeVCMap&
KeModel<T>::getVCMap() {return _impl->getVCMap();}

template<class T>
KeBC<T>&
KeModel<T>::getBC(const int id) {return _impl->getBC(id);}

template<class T>
KeModelOptions<T>&
KeModel<T>::getOptions() {return _impl->getOptions();}


template<class T>
void
KeModel<T>::printBCs()
{
  _impl->printBCs();
}

template<class T>
void
KeModel<T>::advance(const int niter)
{
  _impl->advance(niter);
}

template<class T>
void
KeModel<T>::updateTimek()
{
  _impl->updateTimek();
}

template<class T>
void
KeModel<T>::updateTimee()
{
  _impl->updateTimee();
}


template<class T>
void
KeModel<T>:: getViscosity(const Mesh& mesh)
{
  _impl->getViscosity(mesh);
}
