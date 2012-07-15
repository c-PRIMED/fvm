
#ifndef _COMETMODEL_H_
#define _COMETMODEL_H_

#include <stdio.h>
#include <map>
#include <cmath>
#include <vector>

#ifdef FVM_PARALLEL
  #include <mpi.h>
#endif

#include "Model.h"

#include "Array.h"
#include "Vector.h"

#include "Mesh.h"

#include "Quadrature.h"
#include "DistFunctFields.h"

#include "MacroFields.h"
#include "FlowFields.h"

#include "COMETBC.h"
#include "COMETBoundaryConditions.h"
#include "COMETESBGKDiscretizer.h"

#include "Linearizer.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "Matrix.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"

#include "MatrixOperation.h"
#include "NumType.h"

#include "StressTensor.h"

template<class T>
class COMETModel : public Model
{
 public:
  typedef typename NumTypeTraits<T>:: T_Scalar T_Scalar;
  typedef Array<int> IntArray;
  typedef Array<T> TArray;
  typedef shared_ptr<TArray> TArrptr;
  typedef Array<bool> BArray;
  typedef  Array2D<T> TArray2D;
  typedef Vector<T,3> VectorT3; 
  typedef Array<VectorT3> VectorT3Array;
  typedef shared_ptr<VectorT3Array> VT3Ptr;
  typedef StressTensor<T> StressTensorT6;
  typedef Array<StressTensorT6> StressTensorArray;
  typedef  std::vector<Field*> stdVectorField;
  typedef  DistFunctFields<T> TDistFF;

  typedef Vector<T,5> VectorT5; 
  typedef Array<VectorT5> VectorT5Array;
  typedef Vector<T,6> VectorT6; 
  typedef Array<VectorT6> VectorT6Array;
  typedef Vector<T,10> VectorT10; 
  typedef Array<VectorT10> VectorT10Array;

  typedef shared_ptr<MeshList> MshLstPtr;
  typedef shared_ptr<Mesh> MeshPtr;
  typedef shared_ptr<GeomFields> GeoFldsPtr;
  typedef shared_ptr<StorageSite> SSPtr;
  typedef shared_ptr<CRConnectivity> CRPtr;

  typedef Array<int> BCfaceArray;
  typedef shared_ptr<BCfaceArray> BfacePtr;
  typedef vector<BfacePtr> BCfaceList;
  typedef Array<int> BCcellArray;
  typedef shared_ptr<BCcellArray> BCellPtr;
  typedef vector<BCellPtr> BCcellList;

  typedef Quadrature<T> TQuad;

  typedef std::map<int,COMETBC<T>*> COMETBCMap;
  typedef std::map<int,COMETVC<T>*> COMETVCMap;
  typedef COMETModel<T> TCOMET;
  typedef shared_ptr<TCOMET> TCOMETPtr;

  typedef MultiField::ArrayIndex Index;
  typedef pair<Index,Index> EntryIndex;
  typedef pair<const StorageSite*, const StorageSite*> SSPair;
  
  typedef map<EntryIndex,shared_ptr<Matrix> > MatrixMap;
  typedef map<Index,int> MatrixSizeMap;
  typedef map<const Mesh*,int> SizeMap;
  typedef map<const StorageSite*,StorageSite*> SiteMap;

  typedef map<SSPair,shared_ptr<Array<int> > > MatrixMappersMap;

  typedef map<Index,shared_ptr<StorageSite> > StorageSiteMap;
  typedef map<const StorageSite*,shared_ptr<StorageSite> > GhostStorageSiteMap;

  /**
   * Calculation of macro-parameters density, temperature, components of velocity, pressure
   * by taking moments of distribution function using quadrature points and weights from quadrature.h
   */
  //MacroFields& macroFields;
  
  COMETModel(const MeshList& meshes, const int level, GeomFields& geomFields, MacroFields& macroFields, Quadrature<T>& quad, const int ibm=0,
	     GeomFields* finestGeomFields=NULL, const MeshList* finestMeshes=NULL, MacroFields* finestMacroFields=NULL):
       
    Model(meshes),
    _level(level),
    _geomFields(geomFields),
    _quadrature(quad),
    _macroFields(macroFields),
    _dsfPtr(_meshes,_quadrature,"dsf_"),
    _dsfPtr1(_meshes,_quadrature,"dsf1_"),
    _dsfPtr2(_meshes,_quadrature,"dsf2_"),
    _dsfEqPtr(_meshes,_quadrature,"dsfEq_"),
    _dsfEqPtrES(_meshes,_quadrature,"dsfEqES_"),
    _dsfPtr0(_meshes,_quadrature,"dsf0_"),
    _dsfPtrInj(_meshes,_quadrature,"dsfInj_"),
    _dsfPtrRes(_meshes,_quadrature,"dsfRes_"),
    _dsfPtrFAS(_meshes,_quadrature,"dsfFAS_"),
    _coarseGeomFields("coarse"),
    _initialKmodelNorm(),
    _niters(0),
    _residual(0.0),
    _initialResidual(0.0),
    _ibm(ibm),
    _finestGeomFields(finestGeomFields ? *finestGeomFields : _geomFields),
    _finestMeshes(finestMeshes ? *finestMeshes : _meshes),
    _finestMacroFields(finestMacroFields ? *finestMacroFields : _macroFields)
    {     
     
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
          const StorageSite& faces=mesh.getFaces();
	  const StorageSite& cells=mesh.getCells();
          const int faceCount=faces.getCount();
          const int cellCount=cells.getSelfCount();
          
          BfacePtr BFptr(new BCfaceArray(faceCount));
          BFptr->zero();
          _BFaces.push_back(BFptr);
          
          BCellPtr BCptr(new BCcellArray(cellCount));
          _BCells.push_back(BCptr);
          BCptr->zero();

          BCellPtr ZCptr(new BCcellArray(cellCount));
          _ZCells.push_back(ZCptr);
          ZCptr->zero();
	  if(_level==0)
	  {
	      COMETVC<T> *vc(new COMETVC<T>());
	      vc->vcType = "flow";
	      _vcMap[mesh.getID()] = vc;
	  }
	}

      if(_level==0)
      {
	  SetBoundaryConditions();
	  init();
	  InitializeMacroparameters();
	  initializeMaxwellian();
	  initializeFineMaxwellian();
	  ComputeMacroparameters(); //calculate density,velocity,temperature
	  ComputeFineMacroparameters();
	  ComputeCollisionfrequency(); //calculate viscosity, collisionFrequency
	  initializeMaxwellianEq();    //equilibrium distribution
      }
    }
    
    void init()

    {
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  const Mesh& fMesh = *_finestMeshes[n];

	  const COMETVC<T>& vc = *_vcMap[mesh.getID()];
	  
	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& fCells = fMesh.getCells();	  

	  const int nCells = cells.getCountLevel1();
	  shared_ptr<VectorT3Array> vCell(new VectorT3Array(nCells));
	  
	  VectorT3 initialVelocity;
	  initialVelocity[0] = _options["initialXVelocity"];
	  initialVelocity[1] = _options["initialYVelocity"];
	  initialVelocity[2] = _options["initialZVelocity"];
	  *vCell = initialVelocity;
	  _macroFields.velocity.addArray(cells,vCell);

          shared_ptr<IntArray> fineToCoarseCell(new IntArray(nCells));
          *fineToCoarseCell = -1;
          _geomFields.fineToCoarse.addArray(cells,fineToCoarseCell);

          if((_ibm==1)&&(_level==0))
          {
	      _geomFields.ibType.syncLocal();
              shared_ptr<Array<Vector<int,25> > >finestToCoarseCell(new Array<Vector<int,25> >(fCells.getCountLevel1()));
              Vector<int,25> initialIndex;
              for(int k=0;k<25;k++)
                initialIndex[k]=-1;
              *finestToCoarseCell = initialIndex;
              _finestGeomFields.finestToCoarse.addArray(fCells,finestToCoarseCell);
              Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
              Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);
              for(int c=0;c<nCells;c++)
                FinestToCoarse[c][_level]=c;
          }
	  
	  shared_ptr<TArray> pCell(new TArray(nCells));
	  *pCell = _options["operatingPressure"];
	  _macroFields.pressure.addArray(cells,pCell);
	  
	  shared_ptr<TArray> rhoCell(new TArray(nCells));
	  *rhoCell = 0.;
	  //*rhoCell = vc["density"];
	  _macroFields.density.addArray(cells,rhoCell);
	  
	  shared_ptr<TArray> muCell(new TArray(nCells));
	  *muCell = 0.;
	  //*muCell = vc["viscosity"];
	  _macroFields.viscosity.addArray(cells,muCell);

	  shared_ptr<TArray> tempCell(new TArray(nCells));
	  *tempCell = _options["operatingTemperature"];
	  _macroFields.temperature.addArray(cells,tempCell);

	  shared_ptr<TArray> collFreqCell(new TArray(nCells));
	  *collFreqCell = 0.;
	  //*collFreqCell = vc["viscosity"];
	  _macroFields.collisionFrequency.addArray(cells,collFreqCell);
	  
	  //coeffs for perturbed BGK distribution function
	  shared_ptr<VectorT5Array> coeffCell(new VectorT5Array(nCells));
	  VectorT5 initialCoeff;
	  initialCoeff[0] = 1.0;
	  initialCoeff[1] = 1.0;
	  initialCoeff[2] = 0.0; 
	  initialCoeff[3] = 0.0;
	  initialCoeff[4] = 0.0;
	  *coeffCell = initialCoeff;
	  _macroFields.coeff.addArray(cells,coeffCell);
	  
	  //coeffs for perturbed BGK distribution function
	  shared_ptr<VectorT10Array> coeffgCell(new VectorT10Array(nCells));
	  VectorT10 initialCoeffg;
	  initialCoeffg[0] = 1.0;
	  initialCoeffg[1] = 1.0;
	  initialCoeffg[2] = 0.0; 
	  initialCoeffg[3] = 1.0;
	  initialCoeffg[4] = 0.0; 
	  initialCoeffg[5] = 1.0;
	  initialCoeffg[6] = 0.0;
	  initialCoeffg[7] = 0.0; 
	  initialCoeffg[8] = 0.0;
	  initialCoeffg[9] = 0.0;
	  *coeffgCell = initialCoeffg;
	  _macroFields.coeffg.addArray(cells,coeffgCell);
	  
	  // used for ESBGK equilibrium distribution function
	  shared_ptr<TArray> tempxxCell(new TArray(cells.getCountLevel1()));
	  *tempxxCell = _options["operatingTemperature"]/3;
	  _macroFields.Txx.addArray(cells,tempxxCell);
	  
	  shared_ptr<TArray> tempyyCell(new TArray(cells.getCountLevel1()));
	  *tempyyCell = _options["operatingTemperature"]/3;
	  _macroFields.Tyy.addArray(cells,tempyyCell);
	  
	  shared_ptr<TArray> tempzzCell(new TArray(cells.getCountLevel1()));
	  *tempzzCell = _options["operatingTemperature"]/3;
	  _macroFields.Tzz.addArray(cells,tempzzCell);
	  
	  shared_ptr<TArray> tempxyCell(new TArray(cells.getCountLevel1()));
	  *tempxyCell = 0.0;
	  _macroFields.Txy.addArray(cells,tempxyCell);
	  
	  shared_ptr<TArray> tempyzCell(new TArray(cells.getCountLevel1()));
	  *tempyzCell = 0.0;
	  _macroFields.Tyz.addArray(cells,tempyzCell);
	  
	  shared_ptr<TArray> tempzxCell(new TArray(cells.getCountLevel1()));
	  *tempzxCell = 0.0;
	  _macroFields.Tzx.addArray(cells,tempzxCell);
	  
	//Entropy and Entropy Generation Rate for switching
	  shared_ptr<TArray> EntropyCell(new TArray(cells.getCountLevel1()));
	  *EntropyCell = 0.0;
	  _macroFields.Entropy.addArray(cells,EntropyCell);
	  
	  shared_ptr<TArray> EntropyGenRateCell(new TArray(cells.getCountLevel1()));
	  *EntropyGenRateCell = 0.0;
	  _macroFields.EntropyGenRate.addArray(cells,EntropyGenRateCell);

	  shared_ptr<TArray> EntropyGenRateColl(new TArray(cells.getCountLevel1()));
	  *EntropyGenRateColl = 0.0;
	  _macroFields.EntropyGenRate_Collisional.addArray(cells,EntropyGenRateColl);
	  
	//Pxx,Pyy,Pzz,Pxy,Pxz,Pyz
	shared_ptr<VectorT6Array> stressCell(new VectorT6Array(nCells));
        VectorT6 initialstress;
        initialstress[0] = 1.0;
        initialstress[1] = 1.0;
        initialstress[2] = 1.0; 
	initialstress[3] = 0.0;
	initialstress[4] = 0.0;	
	initialstress[5] = 0.0;
        *stressCell = initialstress;
        _macroFields.Stress.addArray(cells,stressCell);

	//Knq=M300+M120+M102 for Couette with uy
        shared_ptr<TArray> KnqCell(new TArray(cells.getCountLevel1()));
        *KnqCell = 0.0;
        _macroFields.Knq.addArray(cells,KnqCell);


	//higher order moments of distribution function
	/*
	shared_ptr<VectorT3Array> M300Cell(new VectorT3Array(cells.getCount()));
        VectorT3 initialM300;
        initialM300[0] = 0.0;
        initialM300[1] = 0.0;
        initialM300[2] = 0.0; 
        *M300Cell = initialM300;
        _macroFields.M300.addArray(cells,M300Cell);

	shared_ptr<VectorT3Array> M030Cell(new VectorT3Array(cells.getCount()));
        VectorT3 initialM030;
        initialM030[0] = 0.0;
        initialM030[1] = 0.0;
        initialM030[2] = 0.0; 
        *M300Cell = initialM030;
        _macroFields.M030.addArray(cells,M030Cell);
	*/
	//if(MPI::COMM_WORLD.Get_rank()==0)
	//cout<<"array for fields created"<<endl;
	const int numDirections = _quadrature.getDirCount();
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	//FILE * pFile;
	//pFile=fopen("ref_incMEMOSA.txt","w");
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups()){
	  const FaceGroup& fg = *fgPtr; 
	  
	  if((_bcMap[fg.id]->bcType == "SymmetryBC")||(_bcMap[fg.id]->bcType == "RealWallBC")||(_bcMap[fg.id]->bcType == "VelocityInletBC")){
	  
	    const StorageSite& faces = fg.site;
	  		 
	    const Field& areaMagField = _geomFields.areaMag;
	    const TArray& faceAreaMag = dynamic_cast<const TArray &>(areaMagField[faces]);
	    const Field& areaField = _geomFields.area;
	    const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]); 
	  
	      const VectorT3 en = faceArea[0]/faceAreaMag[0];
	      vector<int> tempVec(numDirections);
	      
	      for (int j=0; j<numDirections; j++){
		const T c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
		const T cx_incident = cx[j] - 2.0*c_dot_en*en[0];
		const T cy_incident = cy[j] - 2.0*c_dot_en*en[1];
		const T cz_incident = cz[j] - 2.0*c_dot_en*en[2];           
		int direction_incident=0; 
		T Rdotprod=1e54;
		T dotprod=0.0;
		for (int js=0; js<numDirections; js++){
		  dotprod=pow(cx_incident-cx[js],2)+pow(cy_incident-cy[js],2)+pow(cz_incident-cz[js],2);
		  if (dotprod< Rdotprod){
		    Rdotprod =dotprod;
		    direction_incident=js;}
		}
		tempVec[j] = direction_incident;
		//fprintf(pFile,"%d %d %d \n",fg.id, j,direction_incident);
		
	      }
	      const int fgid=fg.id;
	      _faceReflectionArrayMap[fgid] = tempVec; //add to map
	    
	  }
	}
       
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups()){
	  const FaceGroup& fg = *fgPtr; 
	  if(fg.groupType == "NSinterface"){
	    const StorageSite& Intfaces = fg.site; 
	    
	    shared_ptr<VectorT3Array> InterfaceVelFace(new VectorT3Array(Intfaces.getCount()));
	    InterfaceVelFace ->zero();
	    _macroFields.InterfaceVelocity.addArray(Intfaces,InterfaceVelFace);
	    
	    shared_ptr<StressTensorArray> InterfaceStressFace(new StressTensorArray(Intfaces.getCount()));
	    InterfaceStressFace ->zero();
	    _macroFields.InterfaceStress.addArray(Intfaces,InterfaceStressFace);
	    
	    shared_ptr<TArray> InterfacePressFace(new TArray(Intfaces.getCount()));
	    *InterfacePressFace = _options["operatingPressure"];
	    _macroFields.InterfacePressure.addArray(Intfaces,InterfacePressFace);
	    
	    shared_ptr<TArray> InterfaceDensityFace(new TArray(Intfaces.getCount()));
	    *InterfaceDensityFace =vc["density"];
	    _macroFields.InterfaceDensity.addArray(Intfaces,InterfaceDensityFace);
	    
	    
	  }
	  
	  //fclose(pFile);
	  
	  
	} //end of loop through meshes
      
        BCcellArray& BCArray=*(_BCells[n]);
        BCfaceArray& BCfArray=*(_BFaces[n]);
	BCcellArray& ZCArray=*(_ZCells[n]);
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
            const FaceGroup& fg = *fgPtr;
            if((_bcMap[fg.id]->bcType == "WallBC")||(_bcMap[fg.id]->bcType == "RealWallBC"))
	    {
                const StorageSite& faces = fg.site;
                const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();

                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=2;

                for(int i=0;i<faceCount;i++)
		{
                    int cell1=BfaceCells(i,0);
                    BCArray[cell1]=1;
		}
	    }
            else if(_bcMap[fg.id]->bcType == "VelocityInletBC")
	    {
                const StorageSite& faces = fg.site;
                const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();

                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=3;

                for(int i=0;i<faceCount;i++)
		{
                    int cell1=BfaceCells(i,0);
                    BCArray[cell1]=1;
		}
	    }
            else if(_bcMap[fg.id]->bcType == "ZeroGradBC")
	    {
                const StorageSite& faces = fg.site;
                const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();

                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=4;

                for(int i=0;i<faceCount;i++)
		{
                    int cell1=BfaceCells(i,0);
                    ZCArray[cell1]=1;
		}
	    }
            else if((_bcMap[fg.id]->bcType == "PressureInletBC")||(_bcMap[fg.id]->bcType == "PressureOutletBC"))
	    {
                const StorageSite& faces = fg.site;
                const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();

                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=5;
	    }
            else if(_bcMap[fg.id]->bcType == "SymmetryBC")
	    {
                const StorageSite& faces = fg.site;
                const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();

                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=6;

                for(int i=0;i<faceCount;i++)
		{
                    int cell1=BfaceCells(i,0);
                    BCArray[cell1]=1;
		}
	    }
            else
	    {
                const StorageSite& faces = fg.site;
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();
                         
                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=0;
	    } 
	}
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	{
            const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
	    const int faceCount=faces.getCount();
	    const int offSet=faces.getOffset();
	    
	    for(int i=offSet;i<offSet+faceCount;i++)
	      BCfArray[i]=-1;
	    /*
	    if(MPI::COMM_WORLD.Get_rank()==1)
	      {
		int fC = (mesh.getFaces()).getCount();
		cout<<"level,rank,facecount,iID,offSet,ISize = "<<_level<<" "<<MPI::COMM_WORLD.Get_rank()<<" "<<fC<<" "<<fg.id<<" "<<offSet<<" "<<(offSet+faceCount)<<endl;
	      }
	    */
	}
	
        const StorageSite& faces = mesh.getFaces();
        const int faceCount = faces.getCount();
        const CRConnectivity& faceCells=mesh.getFaceCells(faces);
        const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
        for(int i=0;i<faceCount;i++)
	{
            const int c0 = faceCells(i,0);
            const int c1 = faceCells(i,1);
            if (((ibType[c0] == Mesh::IBTYPE_FLUID) && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
                ((ibType[c1] == Mesh::IBTYPE_FLUID) && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
	    {
                BCfArray[i]=7;
	    }
	}

        int count = 0;
        for(int c=0;c<cells.getSelfCount();c++)
	{
            if(ibType[c] != Mesh::IBTYPE_FLUID)
              count++;
	}
#ifdef FVM_PARALLEL
	MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE, &count, 1, MPI::INT, MPI::SUM);
        if((MPI::COMM_WORLD.Get_rank()==0)&&(_level==0))
          cout<<"number of non-fluid cells in mesh at level "<<_level<<" = "<<count<<endl;
#endif
        
#ifndef FVM_PARALLEL
        if(_level==0)
          cout<<"number of non-fluid cells in mesh at level "<<_level<<" = "<<count<<endl;
#endif
	_niters  =0;
	_initialKmodelNorm = MFRPtr();
	//_initialKmodelvNorm = MFRPtr();
	
	}
    }

    void MakeCoarseModel(TCOMET* finerModel)
    {

      if(_options.AgglomerationMethod=="FaceArea")
      {
	  int maxLevs=finerModel->getOptions().maxLevels;
	  int thisLevel=(finerModel->getLevel())+1;

	  if(thisLevel<maxLevs)  //assumes # of levels will always work for the mesh
	  {
	      MeshList* newMeshesPtr=new MeshList;
	      TQuad* newQuadPtr=new TQuad();
	      MacroFields* newMacroPtr=new MacroFields("coarse");

	      newQuadPtr->CopyQuad(finerModel->getQuadrature());

              MakeCoarseMesh1(finerModel->getMeshList(),
                                           finerModel->getGeomFields(),
                                           *newMeshesPtr);
	      
              _geomFields.fineToCoarse.syncLocal();
	      const Mesh& mesh = *_meshes[0];
	      const StorageSite& cells = mesh.getCells();
	      const int nCells = cells.getCount();
	      Field& FineToCoarseField=(finerModel->getGeomFields()).fineToCoarse;
	      const IntArray& coarseIndex=dynamic_cast<const IntArray&>(FineToCoarseField[cells]);
	      
	      /*	      
	      if(MPI::COMM_WORLD.Get_rank()==1)
		for(int c=0;c<nCells;c++)
		  cout<<" after sync, rank, level, cell no and finetocoarse = "<<MPI::COMM_WORLD.Get_rank()<<" "<<_level<<" "<<c<<" "<<coarseIndex[c]<<endl;
	      */	        
 
	      syncGhostCoarsening(finerModel->getMeshList(),
				  finerModel->getGeomFields(),
				  *newMeshesPtr);
	      const int numMeshes =_meshes.size();
	      for (int n=0; n<numMeshes; n++)
	      {
		  const Mesh& mesh = *_meshes[n];
		  const StorageSite& fineSite = mesh.getCells();
		  StorageSite& coarseSite = *_siteMap[&fineSite];
		  const StorageSite::ScatterMap& fineScatterMap = fineSite.getScatterMap();
		  const StorageSite::ScatterMap& fineScatterMapLevel1 = fineSite.getScatterMapLevel1();
		  StorageSite::ScatterMap& coarseScatterMap = coarseSite.getScatterMap();

		  foreach(const StorageSite::ScatterMap::value_type& pos, fineScatterMap)
		  {
		      const StorageSite& fineOSite = *pos.first;

#ifdef FVM_PARALLEL
		      // the ghost site will not have its corresponding coarse
		      // site created yet so we create it here
		      if (_siteMap.find(&fineOSite) == _siteMap.end())
		      {
			  shared_ptr<StorageSite> ghostSite
			    (new StorageSite(-1));
			  ghostSite->setGatherProcID ( fineOSite.getGatherProcID() );
			  ghostSite->setScatterProcID( fineOSite.getScatterProcID() );
			  ghostSite->setTag( fineOSite.getTag() );
			  StorageSite& coarseOSite = *ghostSite;
			  _siteMap[&fineOSite]=&coarseOSite;
			  _sharedSiteMap[&fineOSite]=ghostSite;
		      }
#endif
   
		      StorageSite& coarseOSite = *_siteMap[&fineOSite];
		      
		      SSPair sskey(&fineSite,&fineOSite);
		      coarseScatterMap[&coarseOSite] = _coarseScatterMaps[sskey];		
		  }

                  foreach(const StorageSite::ScatterMap::value_type& pos, fineScatterMapLevel1)
		  {
                      const StorageSite& fineOSite = *pos.first;
                      SSPair sskey(&fineSite,&fineOSite);
                      if (_coarseScatterMaps.find(sskey) != _coarseScatterMaps.end())
		      {

#ifdef FVM_PARALLEL
                          // the ghost site will not have its corresponding coarse
                          // site created yet so we create it here
                          if (_siteMap.find(&fineOSite) == _siteMap.end())
			  {
                              shared_ptr<StorageSite> ghostSite
                                (new StorageSite(-1));
                              ghostSite->setGatherProcID ( fineOSite.getGatherProcID() );
                              ghostSite->setScatterProcID( fineOSite.getScatterProcID() );
                              ghostSite->setTag( fineOSite.getTag() );
                              StorageSite& coarseOSite = *ghostSite;
                              _siteMap[&fineOSite]=&coarseOSite;
                              _sharedSiteMap[&fineOSite]=ghostSite;
			  }
#endif
			  
                          StorageSite& coarseOSite = *_siteMap[&fineOSite];

                          coarseScatterMap[&coarseOSite] = _coarseScatterMaps[sskey];
		      }
		  }
		  
		  const StorageSite::GatherMap& fineGatherMap = fineSite.getGatherMap();
		  const StorageSite::GatherMap& fineGatherMapLevel1 = fineSite.getGatherMapLevel1();
		  StorageSite::GatherMap& coarseGatherMap = coarseSite.getGatherMap();
		  foreach(const StorageSite::GatherMap::value_type& pos, fineGatherMap)
		  {
		      const StorageSite& fineOSite = *pos.first;
		      StorageSite& coarseOSite = *_siteMap[&fineOSite];
		      SSPair sskey(&fineSite,&fineOSite);

		      coarseGatherMap[&coarseOSite] = _coarseGatherMaps[sskey];
		  }

                  foreach(const StorageSite::GatherMap::value_type& pos, fineGatherMapLevel1)
		  {
                      const StorageSite& fineOSite = *pos.first;
                      SSPair sskey(&fineSite,&fineOSite);
                      if (_coarseGatherMaps.find(sskey) != _coarseGatherMaps.end())
		      {
                          foreach(SiteMap::value_type tempPos, _siteMap)
			  {
                              const StorageSite& tempOSite = *tempPos.first;
                              if(fineOSite.getTag()==tempOSite.getTag())
			      {
                                  //StorageSite& coarseOSite = *_siteMap[&fineOSite];
                                  StorageSite& coarseOSite = *_siteMap[&tempOSite];
                                  coarseGatherMap[&coarseOSite] = _coarseGatherMaps[sskey];
			      }
			  }
		      }
		  }
		  
	      }
	      
	      int newCount= MakeCoarseMesh2(finerModel->getMeshList(),
					    finerModel->getGeomFields(),_coarseGeomFields,
					    *newMeshesPtr);

	      TCOMET* newModelPtr=new COMETModel(*newMeshesPtr,thisLevel,
                                                 _coarseGeomFields,
                                                 *newMacroPtr,*newQuadPtr);
#ifdef FVM_PARALLEL
	      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE, &newCount, 1, MPI::INT, MPI::SUM);
	      if(MPI::COMM_WORLD.Get_rank()==0)
		cout<<"Number of cells in level "<<thisLevel<<"  is "<<newCount<<endl;            
#endif

#ifndef FVM_PARALLEL
	      cout<<"Number of cells in level "<<thisLevel<<"  is "<<newCount<<endl;
#endif

	      newModelPtr->setFinerLevel(finerModel);
	      finerModel->setCoarserLevel(newModelPtr);
	      newModelPtr->getOptions()=finerModel->getOptions();
	      newModelPtr->getBCMap()=finerModel->getBCMap();
	      newModelPtr->getVCMap()=finerModel->getVCMap();

	      newModelPtr->init();
	      newModelPtr->InitializeMacroparameters();
	      newModelPtr->initializeMaxwellian();
	      newModelPtr->initializeCoarseMaxwellian();
	      newModelPtr->ComputeMacroparameters();
	      newModelPtr->ComputeCoarseMacroparameters();
	      newModelPtr->ComputeCollisionfrequency();
	      newModelPtr->initializeMaxwellianEq();

              if(newCount>_options.minCells)
                newModelPtr->MakeCoarseModel(newModelPtr);
              else
                _options.maxLevels=newModelPtr->getLevel();
	  }
      }
      else if(_options.AgglomerationMethod=="AMG")
	throw CException("Have not implemented AMG agglomeration method.");
      else
	throw CException("Unknown agglomeration method.");
    }

    void MakeIBCoarseModel(TCOMET* finerModel, const StorageSite& solidFaces)
    {

      if(_options.AgglomerationMethod=="FaceArea")
        {
          int maxLevs=finerModel->getOptions().maxLevels;
          int thisLevel=(finerModel->getLevel())+1;

          if(thisLevel<maxLevs)  //assumes # of levels will always work for the mesh
            {
              MeshList* newMeshesPtr=new MeshList;
              TQuad* newQuadPtr=new TQuad();
              MacroFields* newMacroPtr=new MacroFields("coarse");

              newQuadPtr->CopyQuad(finerModel->getQuadrature());

              MakeCoarseMesh1(finerModel->getMeshList(),
                              finerModel->getGeomFields(),
                              *newMeshesPtr);
              
              _geomFields.fineToCoarse.syncLocal();
              const Mesh& mesh = *_meshes[0];
              const StorageSite& cells = mesh.getCells();
              const int nCells = cells.getCount();
              Field& FineToCoarseField=(finerModel->getGeomFields()).fineToCoarse;
              const IntArray& coarseIndex=dynamic_cast<const IntArray&>(FineToCoarseField[cells]);
              
              /*              
              if(MPI::COMM_WORLD.Get_rank()==1)
                for(int c=0;c<nCells;c++)
                  cout<<" after sync, rank, level, cell no and finetocoarse = "<<MPI::COMM_WORLD.Get_rank()<<" "<<_level<<" "<<c<<" "<<coarseIndex[c]<<endl;
	      */                
 
              syncGhostCoarsening(finerModel->getMeshList(),
                                  finerModel->getGeomFields(),
                                  *newMeshesPtr);
              const int numMeshes =_meshes.size();
              for (int n=0; n<numMeshes; n++)
                {
                  const Mesh& mesh = *_meshes[n];
                  const StorageSite& fineSite = mesh.getCells();
                  StorageSite& coarseSite = *_siteMap[&fineSite];
                  const StorageSite::ScatterMap& fineScatterMap = fineSite.getScatterMap();
		  const StorageSite::ScatterMap& fineScatterMapLevel1 = fineSite.getScatterMapLevel1();
		  StorageSite::ScatterMap& coarseScatterMap = coarseSite.getScatterMap();

                  foreach(const StorageSite::ScatterMap::value_type& pos, fineScatterMap)
                    {
                      const StorageSite& fineOSite = *pos.first;

#ifdef FVM_PARALLEL
                      // the ghost site will not have its corresponding coarse
                      // site created yet so we create it here
                      if (_siteMap.find(&fineOSite) == _siteMap.end())
                        {
                          shared_ptr<StorageSite> ghostSite
                            (new StorageSite(-1));
                          ghostSite->setGatherProcID ( fineOSite.getGatherProcID() );
                          ghostSite->setScatterProcID( fineOSite.getScatterProcID() );
                          ghostSite->setTag( fineOSite.getTag() );
                          StorageSite& coarseOSite = *ghostSite;
                          _siteMap[&fineOSite]=&coarseOSite;
                          _sharedSiteMap[&fineOSite]=ghostSite;
                        }
#endif
   
                      StorageSite& coarseOSite = *_siteMap[&fineOSite];
                      
                      SSPair sskey(&fineSite,&fineOSite);
                      coarseScatterMap[&coarseOSite] = _coarseScatterMaps[sskey];
                    }
                  
                  foreach(const StorageSite::ScatterMap::value_type& pos, fineScatterMapLevel1)
                  {
                      const StorageSite& fineOSite = *pos.first;
                      SSPair sskey(&fineSite,&fineOSite);
                      if (_coarseScatterMaps.find(sskey) != _coarseScatterMaps.end())
                      {

#ifdef FVM_PARALLEL
                          // the ghost site will not have its corresponding coarse
                          // site created yet so we create it here
                          if (_siteMap.find(&fineOSite) == _siteMap.end())
                          {
                              shared_ptr<StorageSite> ghostSite
                                (new StorageSite(-1));
                              ghostSite->setGatherProcID ( fineOSite.getGatherProcID() );
                              ghostSite->setScatterProcID( fineOSite.getScatterProcID() );
                              ghostSite->setTag( fineOSite.getTag() );
                              StorageSite& coarseOSite = *ghostSite;
                              _siteMap[&fineOSite]=&coarseOSite;
                              _sharedSiteMap[&fineOSite]=ghostSite;
			  }
#endif

                          StorageSite& coarseOSite = *_siteMap[&fineOSite];
                          
                          coarseScatterMap[&coarseOSite] = _coarseScatterMaps[sskey];
		      }
		  }

                  const StorageSite::GatherMap& fineGatherMap = fineSite.getGatherMap();
		  const StorageSite::GatherMap& fineGatherMapLevel1 = fineSite.getGatherMapLevel1();
		  StorageSite::GatherMap& coarseGatherMap = coarseSite.getGatherMap();
                  foreach(const StorageSite::GatherMap::value_type& pos, fineGatherMap)
                    {
                      const StorageSite& fineOSite = *pos.first;
                      StorageSite& coarseOSite = *_siteMap[&fineOSite];
                      SSPair sskey(&fineSite,&fineOSite);

                      coarseGatherMap[&coarseOSite] = _coarseGatherMaps[sskey];
                    }
                  foreach(const StorageSite::GatherMap::value_type& pos, fineGatherMapLevel1)
                  {
                      const StorageSite& fineOSite = *pos.first;
                      SSPair sskey(&fineSite,&fineOSite);
                      if (_coarseGatherMaps.find(sskey) != _coarseGatherMaps.end())
                      {
                          foreach(SiteMap::value_type tempPos, _siteMap)
                          {
                              const StorageSite& tempOSite = *tempPos.first;
                              if(fineOSite.getTag()==tempOSite.getTag())
                              {
                                  //StorageSite& coarseOSite = *_siteMap[&fineOSite];
                                  StorageSite& coarseOSite = *_siteMap[&tempOSite];
                                  coarseGatherMap[&coarseOSite] = _coarseGatherMaps[sskey];
			      }
			  }
		      }
		  } 
                  
                }
              
              int newCount= MakeCoarseMesh2(finerModel->getMeshList(),
                                            finerModel->getGeomFields(),_coarseGeomFields,
                                            *newMeshesPtr);

	      _coarseGeomFields.ibType.syncLocal();

              TCOMET* newModelPtr=new COMETModel(*newMeshesPtr,thisLevel,
                                                 _coarseGeomFields,
                                                 *newMacroPtr,*newQuadPtr,1,&_finestGeomFields,&_finestMeshes,&_finestMacroFields);
#ifdef FVM_PARALLEL
	      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE, &newCount, 1, MPI::INT, MPI::SUM);
              if(MPI::COMM_WORLD.Get_rank()==0)
                cout<<"Number of cells in level "<<thisLevel<<"  is "<<newCount<<endl;
#endif

#ifndef FVM_PARALLEL
              cout<<"Number of cells in level "<<thisLevel<<"  is "<<newCount<<endl;
#endif

              newModelPtr->setFinerLevel(finerModel);
              finerModel->setCoarserLevel(newModelPtr);
              newModelPtr->getOptions()=finerModel->getOptions();
              newModelPtr->getBCMap()=finerModel->getBCMap();
              newModelPtr->getVCMap()=finerModel->getVCMap();

              for (int n=0; n<numMeshes; n++)
	      {
                  const Mesh& mesh = *_meshes[n];
                  const StorageSite& fineIBFaces = mesh.getIBFaces();
		  if(fineIBFaces.getCount()>0)
		  {
		      StorageSite& coarseIBFaces = *_siteMap[&fineIBFaces];
		      for(int dir=0;dir<_quadrature.getDirCount();dir++)
		      {
			  Field& fnd = *_dsfPtr.dsf[dir];
			  const TArray& fIB = dynamic_cast<const TArray&>(fnd[fineIBFaces]);
			  shared_ptr<TArray> cIBV(new TArray(coarseIBFaces.getCount()));
			  cIBV->zero();
			  DistFunctFields<T>& coarserdsf = _coarserLevel->getdsf();
			  Field& cfnd = *coarserdsf.dsf[dir];
			  cfnd.addArray(coarseIBFaces,cIBV);
			  TArray& cIB = dynamic_cast<TArray&>(cfnd[coarseIBFaces]);
			  for(int i=0;i<coarseIBFaces.getCount();i++)
			    cIB[i]=fIB[i];
		      }
		  }

		  shared_ptr<VectorT3Array> coarseSolidVel(new VectorT3Array(solidFaces.getCount()));
		  const VectorT3Array& fineSolidVel = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[solidFaces]);
		  *coarseSolidVel = fineSolidVel;
		  newMacroPtr->velocity.addArray(solidFaces,coarseSolidVel);

                  shared_ptr<TArray> coarseSolidDensity(new TArray(solidFaces.getCount()));
                  const TArray& fineSolidDensity = dynamic_cast<const TArray&>(_macroFields.density[solidFaces]);
                  *coarseSolidDensity = fineSolidDensity;
                  newMacroPtr->density.addArray(solidFaces,coarseSolidDensity);

                  shared_ptr<TArray> coarseSolidTemperature(new TArray(solidFaces.getCount()));
                  const TArray& fineSolidTemperature = dynamic_cast<const TArray&>(_macroFields.temperature[solidFaces]);
                  *coarseSolidTemperature = fineSolidTemperature;
                  newMacroPtr->temperature.addArray(solidFaces,coarseSolidTemperature);
	      }
              
              newModelPtr->init();
              newModelPtr->InitializeMacroparameters();
              newModelPtr->initializeMaxwellian();
              newModelPtr->initializeCoarseMaxwellian();
              newModelPtr->ComputeMacroparameters();
              newModelPtr->ComputeCoarseMacroparameters();
              newModelPtr->ComputeCollisionfrequency();
              newModelPtr->initializeMaxwellianEq();

              if(newCount>_options.minCells)
                newModelPtr->MakeIBCoarseModel(newModelPtr,solidFaces);
              else
                _options.maxLevels=newModelPtr->getLevel();
            }
        }
      else if(_options.AgglomerationMethod=="AMG")
        throw CException("Have not implemented AMG agglomeration method.");
      else
        throw CException("Unknown agglomeration method.");
    }

    void MakeCoarseIndex(const StorageSite& solidFaces)
    {

      const int maxLevs=_options.maxLevels;
      int thisLevel=_level+1;
      
      const Mesh& mesh = *_meshes[0];
      const StorageSite& cells = mesh.getCells();
      const int nCells = cells.getCount();
      const Mesh& fMesh = *_finestMeshes[0];
      const StorageSite& fCells = fMesh.getCells();
      const int nFCells = fCells.getCount();
      Field& FineToCoarseField=_geomFields.fineToCoarse;
      const IntArray& FineToCoarse=dynamic_cast<const IntArray&>(FineToCoarseField[cells]);
      Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
      Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);

      if(_level==0)
        MakeParallel();
      else
      {
	  for(int dir=0;dir<_quadrature.getDirCount();dir++)
	  {
	      Field& fnd = *_dsfPtr.dsf[dir];
	      shared_ptr<TArray> cSV(new TArray(solidFaces.getCount()));
	      cSV->zero();
	      fnd.addArray(solidFaces,cSV);
	  }
      }

      if(thisLevel<maxLevs)  //assumes # of levels will always work for the mesh
      {
          for(int c=0;c<nFCells;c++)
            FinestToCoarse[c][_level+1]=FineToCoarse[FinestToCoarse[c][_level]];
          
          _coarserLevel->MakeCoarseIndex(solidFaces);
      }
    }

    void MakeParallel()
    {

      for(int dir=0;dir<_quadrature.getDirCount();dir++)
	{
          Field& fnd = *_dsfPtr.dsf[dir];
          fnd.syncLocal();
	}
    }
    
  void InitializeMacroparameters()
  {  const int numMeshes =_meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount(); 
	//	const double pi(3.14159);
	const double pi=_options.pi;
	TArray& Entropy = dynamic_cast<TArray&>(_macroFields.Entropy[cells]);  
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);  
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[cells]);
	
	VectorT5Array& coeff = dynamic_cast<VectorT5Array&>(_macroFields.coeff[cells]);
	VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
	TArray& Txx = dynamic_cast<TArray&>(_macroFields.Txx[cells]);
	TArray& Tyy = dynamic_cast<TArray&>(_macroFields.Tyy[cells]);
	TArray& Tzz = dynamic_cast<TArray&>(_macroFields.Tzz[cells]);
	TArray& Txy = dynamic_cast<TArray&>(_macroFields.Txy[cells]);
	TArray& Tyz = dynamic_cast<TArray&>(_macroFields.Tyz[cells]);
	TArray& Tzx = dynamic_cast<TArray&>(_macroFields.Tzx[cells]);
	//if ( MPI::COMM_WORLD.Get_rank() == 0 ) {cout << "ncells="<<nCells<<endl;}
	
	TArray& Knq = dynamic_cast<TArray&>(_macroFields.Knq[cells]);
	//initialize density,velocity  
	for(int c=0; c<nCells;c++)
	  {
	    density[c] =1.0;
	    v[c][0]=0.0;
	    v[c][1]=0.0;
	    v[c][2]=0.0;
	    temperature[c]=1.0;
	    pressure[c]=temperature[c]*density[c];
	    
	    //BGK
	      coeff[c][0]=density[c]/pow((pi*temperature[c]),1.5);
	      coeff[c][1]=1/temperature[c];
	      coeff[c][2]=0.0;coeff[c][3]=0.0;coeff[c][4]=0.0;
	    
	    Entropy[c]=0.0;
            if(_options.fgamma ==2){
	    //ESBGK
	    coeffg[c][0]=coeff[c][0];
	    coeffg[c][1]=coeff[c][1];
	    coeffg[c][2]=coeff[c][2];
	    coeffg[c][3]=coeff[c][1];
	    coeffg[c][4]=coeff[c][3];
	    coeffg[c][5]=coeff[c][1];
	    coeffg[c][6]=coeff[c][4];
	    coeffg[c][7]=0.0;
	    coeffg[c][8]=0.0;
	    coeffg[c][9]=0.0;	
	    }    
	    
	    Txx[c]=0.5;
	    Tyy[c]=0.5;
	    Tzz[c]=0.5;
	    Txy[c]=0.0;
	    Tyz[c]=0.0;
	    Tzx[c]=0.0;
	    
	    Knq[c]=0.0;
	  }	
      }
  }
  void InitializeFgammaCoefficients()
  {
    const int numMeshes =_meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount(); 
	//	const double pi(3.14159);
	const double pi=_options.pi;
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);  
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	VectorT5Array& coeff = dynamic_cast<VectorT5Array&>(_macroFields.coeff[cells]);
	VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
	TArray& Txx = dynamic_cast<TArray&>(_macroFields.Txx[cells]);
	TArray& Tyy = dynamic_cast<TArray&>(_macroFields.Tyy[cells]);
	TArray& Tzz = dynamic_cast<TArray&>(_macroFields.Tzz[cells]);
	TArray& Txy = dynamic_cast<TArray&>(_macroFields.Txy[cells]);
	TArray& Tyz = dynamic_cast<TArray&>(_macroFields.Tyz[cells]);
	TArray& Tzx = dynamic_cast<TArray&>(_macroFields.Tzx[cells]);
	
	for(int c=0; c<nCells;c++)
	  {
	    
	      //BGK
	      coeff[c][0]=density[c]/pow((pi*temperature[c]),1.5);
	      coeff[c][1]=1/temperature[c];
	      coeff[c][2]=0.0;coeff[c][3]=0.0;coeff[c][4]=0.0;
	   
	    if(_options.fgamma ==2){
	      //ESBGK
	      coeffg[c][0]=coeff[c][0];
	      coeffg[c][1]=coeff[c][1];
	      coeffg[c][2]=coeff[c][2];
	      coeffg[c][3]=coeff[c][1];
	      coeffg[c][4]=coeff[c][3];
	      coeffg[c][5]=coeff[c][1];
	      coeffg[c][6]=coeff[c][4];
	      coeffg[c][7]=0.0;
	      coeffg[c][8]=0.0;
	      coeffg[c][9]=0.0;	

	      Txx[c]=0.5*temperature[c];
	      Tyy[c]=0.5*temperature[c];
	      Tzz[c]=0.5*temperature[c];
	      Txy[c]=0.0;
	      Tyz[c]=0.0;
	      Tzx[c]=0.0;
	    }
	  }
      }
  }

  void ComputeMacroparameters() 
  {  
    //FILE * pFile;
    //pFile = fopen("distfun_mf.txt","w");
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();  //
	
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[cells]);
	const int N123 = _quadrature.getDirCount(); 
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	
	VectorT6Array& stress = dynamic_cast<VectorT6Array&>(_macroFields.Stress[cells]);
	
	//initialize density,velocity,temperature to zero    
	for(int c=0; c<nCells;c++)
	  {
	    density[c]=0.0;
	    v[c][0]=0.0;
	    v[c][1]=0.0;
	    v[c][2]=0.0;
	    temperature[c]=0.0;   
	    stress[c][0]=0.0;stress[c][1]=0.0;stress[c][2]=0.0;
	    stress[c][3]=0.0;stress[c][4]=0.0;stress[c][5]=0.0;

	  }	


	for(int j=0;j<N123;j++){
	  
	  Field& fnd = *_dsfPtr.dsf[j];
	  const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
	  
	  //fprintf(pFile,"%d %12.6f %E %E %E %E \n",j,dcxyz[j],cx[j],cy[j],f[80],density[80]+dcxyz[j]*f[80]);
	  
	  for(int c=0; c<nCells;c++){
	    density[c] = density[c]+wts[j]*f[c];
	    v[c][0]= v[c][0]+(cx[j]*f[c])*wts[j];
	    v[c][1]= v[c][1]+(cy[j]*f[c])*wts[j];
	    v[c][2]= v[c][2]+(cz[j]*f[c])*wts[j];
	    temperature[c]= temperature[c]+(pow(cx[j],2.0)+pow(cy[j],2.0)
					   +pow(cz[j],2.0))*f[c]*wts[j];
	   
	  }
	  
	}
	


	for(int c=0; c<nCells;c++){
	  v[c][0]=v[c][0]/density[c];
	  v[c][1]=v[c][1]/density[c];
	  v[c][2]=v[c][2]/density[c];
	  temperature[c]=temperature[c]-(pow(v[c][0],2.0)
					 +pow(v[c][1],2.0)
					 +pow(v[c][2],2.0))*density[c];
	  temperature[c]=temperature[c]/(1.5*density[c]);  
	  pressure[c]=density[c]*temperature[c];
	}
	
	//Find Pxx,Pyy,Pzz,Pxy,Pyz,Pzx, etc in field	
	
	for(int j=0;j<N123;j++){	  
	  Field& fnd = *_dsfPtr.dsf[j];
	  const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);	  
	  for(int c=0; c<nCells;c++){
	    stress[c][0] +=pow((cx[j]-v[c][0]),2.0)*f[c]*wts[j];
	    stress[c][1] +=pow((cy[j]-v[c][1]),2.0)*f[c]*wts[j];
	    stress[c][2] +=pow((cz[j]-v[c][2]),2.0)*f[c]*wts[j];
	    stress[c][3] +=(cx[j]-v[c][0])*(cy[j]-v[c][1])*f[c]*wts[j];
	    stress[c][4] +=(cy[j]-v[c][1])*(cz[j]-v[c][2])*f[c]*wts[j];
	    stress[c][5] +=(cz[j]-v[c][2])*(cx[j]-v[c][0])*f[c]*wts[j];
	    
	  }}
	
	
      }// end of loop over nmeshes
    //fclose(pFile);
  }

  void ComputeCOMETMacroparameters() 
  {  
    //FILE * pFile;
    //pFile = fopen("distfun_mf.txt","w");
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {

        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const int nCells = cells.getCount();  //


        TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);
        TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);
        const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
        TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[cells]);
        const int N123 = _quadrature.getDirCount(); 

        const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
        const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
        const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
        const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

        VectorT6Array& stress = dynamic_cast<VectorT6Array&>(_macroFields.Stress[cells]);

        //initialize density,velocity,temperature to zero    
        for(int c=0; c<nCells;c++)
          {
            density[c]=0.0;
            temperature[c]=0.0;   
            stress[c][0]=0.0;stress[c][1]=0.0;stress[c][2]=0.0;
            stress[c][3]=0.0;stress[c][4]=0.0;stress[c][5]=0.0;

          }


        for(int j=0;j<N123;j++){
          
          Field& fnd = *_dsfPtr.dsf[j];
          const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
          
          //fprintf(pFile,"%d %12.6f %E %E %E %E \n",j,dcxyz[j],cx[j],cy[j],f[80],density[80]+dcxyz[j]*f[80]);
          
          for(int c=0; c<nCells;c++){
            density[c] = density[c]+wts[j]*f[c];
            temperature[c]= temperature[c]+(pow(cx[j],2.0)+pow(cy[j],2.0)
					    +pow(cz[j],2.0))*f[c]*wts[j];
           
          }
          
        }



        for(int c=0; c<nCells;c++){
          temperature[c]=temperature[c]-(pow(v[c][0],2.0)
                                         +pow(v[c][1],2.0)
                                         +pow(v[c][2],2.0))*density[c];
          temperature[c]=temperature[c]/(1.5*density[c]);  
          pressure[c]=density[c]*temperature[c];
        }

        //Find Pxx,Pyy,Pzz,Pxy,Pyz,Pzx, etc in field

        for(int j=0;j<N123;j++){          
          Field& fnd = *_dsfPtr.dsf[j];
          const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);      
          for(int c=0; c<nCells;c++){
            stress[c][0] +=pow((cx[j]-v[c][0]),2.0)*f[c]*wts[j];
            stress[c][1] +=pow((cy[j]-v[c][1]),2.0)*f[c]*wts[j];
            stress[c][2] +=pow((cz[j]-v[c][2]),2.0)*f[c]*wts[j];
            stress[c][3] +=(cx[j]-v[c][0])*(cy[j]-v[c][1])*f[c]*wts[j];
            stress[c][4] +=(cy[j]-v[c][1])*(cz[j]-v[c][2])*f[c]*wts[j];
            stress[c][5] +=(cz[j]-v[c][2])*(cx[j]-v[c][0])*f[c]*wts[j];
            
          }}


      }// end of loop over nmeshes
    //fclose(pFile);
  }

  void ComputeFineMacroparameters() 
  {  
    const int numMeshes = _meshes.size();
    const T zero(0.0);
    for (int n=0; n<numMeshes; n++)
      {

        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const int nCells = cells.getCount();  //

	VectorT3 zeroVelocity;
	zeroVelocity[0] = zero;
	zeroVelocity[1] = zero;
	zeroVelocity[2] = zero;

	shared_ptr<VectorT3Array> vRCell(new VectorT3Array(nCells));
	*vRCell = zeroVelocity;
	_macroFields.velocityResidual.addArray(cells,vRCell);

	/*
        VectorT3Array& vR = dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[cells]);
	
	for(int c=0; c<nCells;c++)
          {
            vR[c][0]=0.0;
            vR[c][1]=0.0;
            vR[c][2]=0.0;
          }
	*/

      }// end of loop over nmeshes
    //fclose(pFile);
  }

  void ComputeCoarseMacroparameters()
  {
    const int numMeshes = _meshes.size();
    const T zero(0.0);
    for (int n=0; n<numMeshes; n++)
      {

        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const int nCells = cells.getCount();  //

	VectorT3 zeroVelocity;
        zeroVelocity[0] = zero;
	zeroVelocity[1] = zero;
	zeroVelocity[2] = zero;

	shared_ptr<VectorT3Array> vRCell(new VectorT3Array(nCells));
        *vRCell = zeroVelocity;
        _macroFields.velocityResidual.addArray(cells,vRCell);

        shared_ptr<VectorT3Array> vICell(new VectorT3Array(nCells));
        *vICell = zeroVelocity;
        _macroFields.velocityInjected.addArray(cells,vICell);

        shared_ptr<VectorT3Array> vFCell(new VectorT3Array(nCells));
        *vFCell = zeroVelocity;
        _macroFields.velocityFASCorrection.addArray(cells,vFCell);

	/*
        VectorT3Array& vR = dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[cells]);
        VectorT3Array& vI = dynamic_cast<VectorT3Array&>(_macroFields.velocityInjected[cells]);
        VectorT3Array& vF = dynamic_cast<VectorT3Array&>(_macroFields.velocityFASCorrection[cells]);

        for(int c=0; c<nCells;c++)
          {
            vR[c][0]=0.0;
            vR[c][1]=0.0;
            vR[c][2]=0.0;

            vI[c][0]=0.0;
            vI[c][1]=0.0;
            vI[c][2]=0.0;

            vF[c][0]=0.0;
            vF[c][1]=0.0;
            vF[c][2]=0.0;
          }
	*/

      }// end of loop over nmeshes
    //fclose(pFile);
  }
  
  void ComputeMacroparametersESBGK() 
  {  
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {	
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);

	TArray& Txx = dynamic_cast<TArray&>(_macroFields.Txx[cells]);
	TArray& Tyy = dynamic_cast<TArray&>(_macroFields.Tyy[cells]);
	TArray& Tzz = dynamic_cast<TArray&>(_macroFields.Tzz[cells]);
	TArray& Txy = dynamic_cast<TArray&>(_macroFields.Txy[cells]);
	TArray& Tyz = dynamic_cast<TArray&>(_macroFields.Tyz[cells]);
	TArray& Tzx = dynamic_cast<TArray&>(_macroFields.Tzx[cells]);
	
	const int N123 = _quadrature.getDirCount(); 
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	
	const double Pr=_options.Prandtl;
	//cout <<"Prandlt" <<Pr<<endl;
	
	//initialize density,velocity,temperature to zero    
	for(int c=0; c<nCells;c++)
	  {
	    Txx[c]=0.0;
	    Tyy[c]=0.0;
	    Tzz[c]=0.0;
	    Txy[c]=0.0;
	    Tyz[c]=0.0;
	    Tzx[c]=0.0;
	  }
	for(int j=0;j<N123;j++){
	  
	  Field& fnd = *_dsfPtr.dsf[j];
	  const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
	  Field& fndEq = *_dsfEqPtr.dsf[j];
	  const TArray& fgam = dynamic_cast<const TArray&>(fndEq[cells]);	  
	  for(int c=0; c<nCells;c++){
	    Txx[c]=Txx[c]+pow(cx[j]-v[c][0],2)*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];
	    Tyy[c]=Tyy[c]+pow(cy[j]-v[c][1],2)*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j] ;
	    Tzz[c]=Tzz[c]+pow(cz[j]-v[c][2],2)*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];
	    //Txy[c]=Txy[c]+(cx[j])*(cy[j])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];
	    //Tyz[c]=Tyz[c]+(cy[j])*(cz[j])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];            //Tzx[c]=Tzx[c]+(cz[j])*(cx[j])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];

 
	    Txy[c]=Txy[c]+(cx[j]-v[c][0])*(cy[j]-v[c][1])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];
	    Tyz[c]=Tyz[c]+(cy[j]-v[c][1])*(cz[j]-v[c][2])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j]; 
	    Tzx[c]=Tzx[c]+(cz[j]-v[c][2])*(cx[j]-v[c][0])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j]; 
	  }
	}
		
	for(int c=0; c<nCells;c++){
	  Txx[c]=Txx[c]/density[c];
	  Tyy[c]=Tyy[c]/density[c];
	  Tzz[c]=Tzz[c]/density[c];
	  Txy[c]=Txy[c]/density[c];
	  Tyz[c]=Tyz[c]/density[c];
	  Tzx[c]=Tzx[c]/density[c];
	}

     }
  }
  
  

  /*
   * Collision frequency
   *
   *
   */

  void ComputeCollisionfrequency()  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	
	const T rho_init=_options["rho_init"]; 
	const T T_init= _options["T_init"]; 
	const T mu_w= _options["mu_w"];
	const T Tmuref= _options["Tmuref"];
	const T muref= _options["muref"];
	const T R=8314.0/_options["molecularWeight"];
	const T nondim_length=_options["nonDimLt"];

	const T mu0=rho_init*R* T_init*nondim_length/pow(2*R* T_init,0.5);  
	
	  
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);
	TArray& viscosity = dynamic_cast<TArray&>(_macroFields.viscosity[cells]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);

	TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
	for(int c=0; c<nCells;c++)
	  {
	    viscosity[c]= muref*pow(temperature[c]*T_init/ Tmuref,mu_w); // viscosity power law
	    collisionFrequency[c]=density[c]*temperature[c]/viscosity[c]*mu0;
	  }
	
	if(_options.fgamma==2){
	  for(int c=0; c<nCells;c++)
	    collisionFrequency[c]=_options.Prandtl*collisionFrequency[c];
	}
	
      }
  }
  
  void MomentHierarchy()  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const int Knq_dir=_options.Knq_direction; 
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	TArray& Knq = dynamic_cast<TArray&>(_macroFields.Knq[cells]);
	const int num_directions = _quadrature.getDirCount(); 
	if (Knq_dir ==0){
	  for(int j=0;j<num_directions;j++){
	    Field& fnd = *_dsfPtr.dsf[j];
	    const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
	    for(int c=0; c<nCells;c++){
	      Knq[c]=Knq[c]+0.5*f[c]*wts[j]*(pow(cx[j]-v[c][0],3.0)+(cx[j]-v[c][0])*pow(cy[j]-v[c][1],2.0)+(cx[j]-v[c][0])*pow(cz[j]-v[c][2],2.0));
	  }
	  }}
	else if(Knq_dir ==1){
	  for(int j=0;j<num_directions;j++){
	    Field& fnd = *_dsfPtr.dsf[j];
	    const TArray& f = dynamic_cast<const TArray&>(fnd[cells]); 
	    for(int c=0; c<nCells;c++){
	      Knq[c]=Knq[c]+0.5*f[c]*wts[j]*(pow(cy[j]-v[c][1],3.0)+(cy[j]-v[c][1])*pow(cx[j]-v[c][0],2.0)+(cy[j]-v[c][1])*pow(cz[j]-v[c][2],2.0));
	    }
	  }}
	
	else if(Knq_dir ==2){
	  for(int j=0;j<num_directions;j++){
	    Field& fnd = *_dsfPtr.dsf[j];
	    const TArray& f = dynamic_cast<const TArray&>(fnd[cells]); 
	    for(int c=0; c<nCells;c++){
	      Knq[c]=Knq[c]+0.5*f[c]*wts[j]*(pow(cz[j]-v[c][2],3.0)+(cz[j]-v[c][2])*pow(cx[j]-v[c][0],2.0)+(cz[j]-v[c][2])*pow(cy[j]-v[c][1],2.0));
	    }
	  }}
      }
  }
  void EntropyGeneration()  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const T rho_init=_options["rho_init"]; 
	const T T_init= _options["T_init"]; 
	const T molwt=_options["molecularWeight"]*1E-26/6.023;
	const T R=8314.0/_options["molecularWeight"];
	const T u_init=pow(2.0*R*T_init,0.5);
	const T Planck=_options.Planck;
	const T h3bm4u3=pow(Planck,3)/ pow(molwt,4)*rho_init/pow(u_init,3);
	//cout << "h3bm4u3  " << h3bm4u3 <<endl;	
	//cout <<" u_init "<<u_init<<" rho_init "<<rho_init<<endl;
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);	
	TArray& Entropy = dynamic_cast<TArray&>(_macroFields.Entropy[cells]);
	TArray& EntropyGenRate_Collisional = dynamic_cast<TArray&>(_macroFields.EntropyGenRate_Collisional[cells]);
	TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
	for(int c=0; c<nCells;c++){ 
	  Entropy[c]=0.0;EntropyGenRate_Collisional[c]=0.0;
	}
	const int num_directions = _quadrature.getDirCount(); 
	if (_options.fgamma ==2){
	  for(int j=0;j<num_directions;j++){
	    Field& fnd = *_dsfPtr.dsf[j];
	    Field& feqES = *_dsfEqPtrES.dsf[j]; //for fgamma_2
	    const TArray& f = dynamic_cast<const TArray&>(fnd[cells]); 
	    const TArray& fgam = dynamic_cast<const TArray&>(feqES[cells]);
	    for(int c=0; c<nCells;c++){
	      Entropy[c]=Entropy[c]+f[c]*wts[j]*(1-log(h3bm4u3*f[c]));
	      EntropyGenRate_Collisional[c]+= (f[c]-fgam[c])*collisionFrequency[c]*log(h3bm4u3*f[c])*wts[j];
	    }
	  }
	}
	else{
	  for(int j=0;j<num_directions;j++){
	    Field& fnd = *_dsfPtr.dsf[j];
	    Field& feq = *_dsfEqPtr.dsf[j];
	    const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
	    const TArray& fgam = dynamic_cast<const TArray&>(feq[cells]);
	    for(int c=0; c<nCells;c++){
	      Entropy[c]=Entropy[c]+f[c]*wts[j]*(1-log(h3bm4u3*f[c]));
	      EntropyGenRate_Collisional[c]+=(f[c]-fgam[c])*collisionFrequency[c]*(log(h3bm4u3*f[c]))*wts[j];
	    }
	  }
	}
	
	
      }
  }
  
  void initializeMaxwellianEq()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	const VectorT5Array& coeff = dynamic_cast<VectorT5Array&>(_macroFields.coeff[cells]);
	const VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 

	for(int j=0;j< numFields;j++){
	  Field& fndEq = *_dsfEqPtr.dsf[j];
	  TArray& fEq = dynamic_cast< TArray&>(fndEq[cells]);
	  for(int c=0; c<nCells;c++){
	    fEq[c]=coeff[c][0]*exp(-coeff[c][1]*(pow(cx[j]-v[c][0],2)+pow(cy[j]-v[c][1],2)
						 +pow(cz[j]-v[c][2],2))+coeff[c][2]*(cx[j]-v[c][0])
				   +coeff[c][3]*(cy[j]-v[c][1])+coeff[c][4]*(cz[j]-v[c][2]));
	  } 
	}
	
	if(_options.fgamma==2){
	  
	  for(int j=0;j< numFields;j++){
	    Field& fndEqES = *_dsfEqPtrES.dsf[j];
	    TArray& fEqES = dynamic_cast< TArray&>(fndEqES[cells]);
	    for(int c=0; c<nCells;c++){ 
	      T Cc1=(cx[j]-v[c][0]);
	      T Cc2=(cy[j]-v[c][1]);
	      T Cc3=(cz[j]-v[c][2]);
	      fEqES[c]=coeffg[c][0]*exp(-coeffg[c][1]*pow(Cc1,2)+coeffg[c][2]*Cc1
					-coeffg[c][3]*pow(Cc2,2)+coeffg[c][4]*Cc2
					-coeffg[c][5]*pow(Cc3,2)+coeffg[c][6]*Cc3
					+coeffg[c][7]*cx[j]*cy[j]+coeffg[c][8]*cy[j]*cz[j]
					+coeffg[c][9]*cz[j]*cx[j]);
	    }
	  }
	}
	
	
	
      }
  }
  
  void NewtonsMethodBGK(const int ktrial)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {		
	//cout << " NewtonsMethod" <<endl;
	const T tolx=_options["ToleranceX"];
	const T tolf=_options["ToleranceF"];
	const int sizeC=5;
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
    	
	VectorT5Array& coeff = dynamic_cast<VectorT5Array&>(_macroFields.coeff[cells]);
	
	for(int c=0; c<nCells;c++){
	  
	  for (int trial=0;trial<ktrial;trial ++){
	    SquareMatrix<T,sizeC> fjac(0);
	    SquareMatrix<T,sizeC> fjacinv(0);
	    VectorT5 fvec;
	   
	    
	    fvec[0]=density[c];
	    fvec[1]=density[c]*v[c][0];
	    fvec[2]=density[c]*v[c][1];
	    fvec[3]=density[c]*v[c][2];
	    fvec[4]=1.5*density[c]*temperature[c]+density[c]*(pow(v[c][0],2)+pow(v[c][1],2)+pow(v[c][2],2.0));
	    
	  
	    setJacobianBGK(fjac,fvec,coeff[c],v[c],c);
	    
	    //solve using GE or inverse
	    T errf=0.;
	    for (int row=0;row<sizeC;row++){errf+=fabs(fvec[row]);}
	    
	    if(errf <= tolf)
	      break;
	    VectorT5 pvec;
	    for (int row=0;row<sizeC;row++){pvec[row]=-fvec[row];}//rhs
	  

	    //solve Ax=b for x 
	    //p=GE_elim(fjac,p,3);
	    VectorT5 xvec;
	    fjacinv=inverseGauss(fjac,sizeC);
	    
	   
	    
	    for (int row=0;row<sizeC;row++){ 
	      xvec[row]=0.0;
	      for (int col=0;col<sizeC;col++){
		xvec[row]+=fjacinv(row,col)*pvec[col];}
	    }
	    //check for convergence, update
	    
	    T errx=0.;
	    for (int row=0;row<sizeC;row++){
	      errx +=fabs(xvec[row]);
	      coeff[c][row]+= xvec[row];
	    }
	    
	  
	    if(errx <= tolx)
	      break;
	    
	  }
	  
	  
	}
      }
    
  }
  
  void EquilibriumDistributionBGK()
  {	
    const int  ktrial=_options.NewtonsMethod_ktrial;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	//const double pi=_options.pi;
	//const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	//const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
	
	//initialize coeff
	VectorT5Array& coeff = dynamic_cast<VectorT5Array&>(_macroFields.coeff[cells]);
	
	/*
	for(int c=0; c<nCells;c++){
	  coeff[c][0]=density[c]/pow((pi*temperature[c]),1.5);
	  coeff[c][1]=1/temperature[c];
	  coeff[c][2]=0.0;
	  coeff[c][3]=0.0;
	  coeff[c][4]=0.0;	    
	}
	*/
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	//call Newtons Method
	
	NewtonsMethodBGK(ktrial);

	//calculate perturbed maxwellian for BGK
	for(int j=0;j< numFields;j++){
	  Field& fndEq = *_dsfEqPtr.dsf[j];
	  TArray& fEq = dynamic_cast< TArray&>(fndEq[cells]);
	  for(int c=0; c<nCells;c++){
	    fEq[c]=coeff[c][0]*exp(-coeff[c][1]*(pow(cx[j]-v[c][0],2)+pow(cy[j]-v[c][1],2)
	    +pow(cz[j]-v[c][2],2))+coeff[c][2]*(cx[j]-v[c][0])
	      			    +coeff[c][3]*(cy[j]-v[c][1])+coeff[c][4]*(cz[j]-v[c][2]));
	  } 
	  
	}
      }
  }
  
  
  void setJacobianBGK(SquareMatrix<T,5>& fjac, VectorT5& fvec, const VectorT5& xn,const VectorT3& v,const int c)
  { 
    
   
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const TArray& wts = dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);	
    
    const TArray2D& malphaBGK = dynamic_cast<const TArray2D&>(*_quadrature.malphaBGKPtr);	
    const int numFields= _quadrature.getDirCount(); 
    VectorT5 mexp;
    
    for(int j=0;j< numFields;j++){   
      T Cconst=pow(cx[j]-v[0],2.0)+pow(cy[j]-v[1],2.0)+pow(cz[j]-v[2],2.0);
      T Econst=xn[0]*exp(-xn[1]*Cconst+xn[2]*(cx[j]-v[0])+xn[3]*(cy[j]-v[1])+xn[4]*(cz[j]-v[2]))*wts[j];
      
      for (int row=0;row<5;row++){
	fvec[row]+= -Econst*malphaBGK(j,row); //smm
	
	//fvec[row]=tvec[row]+fvec[row];               //mma  
      }
     
      mexp[0]=-Econst/xn[0];
      mexp[1]=Econst*Cconst;
      mexp[2]=-Econst*(cx[j]-v[0]);
      mexp[3]=-Econst*(cy[j]-v[1]);
      mexp[4]=-Econst*(cz[j]-v[2]);
      for (int row=0;row<5;row++){
	for (int col=0;col<5;col++){
	  fjac(row,col)+=malphaBGK(j,row)*mexp[col];  //new
	}
      }
      
      
    }
    
    
  
  }
  
  void NewtonsMethodESBGK(const int ktrial)
  {
    // cout<< "Inside Newtons Method" <<endl;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const T tolx=_options["ToleranceX"];
	const T tolf=_options["ToleranceF"];
	const int sizeC=10;
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	
	const TArray& Txx = dynamic_cast<const TArray&>(_macroFields.Txx[cells]);
	const TArray& Tyy = dynamic_cast<const TArray&>(_macroFields.Tyy[cells]);
	const TArray& Tzz = dynamic_cast<const TArray&>(_macroFields.Tzz[cells]);
	const TArray& Txy = dynamic_cast<const TArray&>(_macroFields.Txy[cells]);
	const TArray& Tyz = dynamic_cast<const TArray&>(_macroFields.Tyz[cells]);
	const TArray& Tzx = dynamic_cast<const TArray&>(_macroFields.Tzx[cells]);
	
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
    	
	VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
	
	for(int c=0; c<nCells;c++){
	  // {cout <<"NM:ES" <<c <<endl;}
	  for (int trial=0;trial<ktrial;trial ++){
	    SquareMatrix<T,sizeC> fjac(0);
	    SquareMatrix<T,sizeC> fjacinv(0);
	    Vector<T,sizeC> fvec;
	    // if (c==_options.printCellNumber){cout <<"trial" <<trial <<endl;}
	    
	    fvec[0]=density[c];
	    fvec[1]=density[c]*v[c][0];
	    fvec[2]=density[c]*v[c][1];
	    fvec[3]=density[c]*v[c][2];
	    fvec[4]=density[c]*(pow(v[c][0],2)+Txx[c]); 
	    fvec[5]=density[c]*(pow(v[c][1],2)+Tyy[c]);
	    fvec[6]=density[c]*(pow(v[c][2],2)+Tzz[c]);
	    fvec[7]=density[c]*(v[c][0]*v[c][1]+Txy[c]);
	    fvec[8]=density[c]*(v[c][1]*v[c][2]+Tyz[c]);
	    fvec[9]=density[c]*(v[c][2]*v[c][0]+Tzx[c]);
	    
	    //calculate Jacobian
	    setJacobianESBGK(fjac,fvec,coeffg[c],v[c],c);
	    
	    
	    //solve using GaussElimination
	    T errf=0.; //and Jacobian matrix in fjac.
	    for (int row=0;row<sizeC;row++){errf+=fabs(fvec[row]);}
	    if(errf <= tolf)
	      break;
	    Vector<T,sizeC> pvec;
	    for (int row=0;row<sizeC;row++){pvec[row]=-fvec[row];}
	    
	    
	    //solve Ax=b for x 
	    //p=GE_elim(fjac,p,3);
	    Vector<T,sizeC> xvec; 
	    fjacinv=inverseGauss(fjac,sizeC);
	    for (int row=0;row<sizeC;row++){ 
	      xvec[row]=0.0;
	      for (int col=0;col<sizeC;col++){
		xvec[row]+=fjacinv(row,col)*pvec[col];
	      
	      }
	    }
	    /*
	    if (c==_options.printCellNumber){
	      cout << " cg0 "<<coeffg[c][0]<<" cg1 "<<coeffg[c][1]<<"  cg2 "<<coeffg[c][2] << endl;
	      cout <<" cg3 " <<coeffg[c][3]<< " cg4 "<<coeffg[c][4]<<" cg5 "<<coeffg[c][5]<<"  cg6 "<<coeffg[c][6] << endl;
	      cout <<" cg7 " <<coeffg[c][7]<< " cg8 "<<coeffg[c][8]<<" cg9 "<<coeffg[c][9]<<endl;
	    
	          
	    //cout << " fvec-ESBGK " << fvec[4] <<fvec[5]<<fvec[6]<<fvec[7]<<fvec[8]<<fvec[9] <<endl;
	    FILE * pFile;
	    pFile = fopen("fvecfjac.dat","wa");  
	    //fprintf(pFile,"%s %d \n","trial",trial);
	    for (int mat_col=0;mat_col<sizeC;mat_col++){fprintf(pFile,"%12.4E",fvec[mat_col]);} 
	    fprintf(pFile,"\n");
	    for (int mat_row=0;mat_row<sizeC;mat_row++){
	      for (int mat_col=0;mat_col<sizeC;mat_col++){
		  fprintf(pFile,"%12.4E",fjac(mat_row,mat_col));}
	      fprintf(pFile,"\n");} 
	    // fprintf(pFile,"done \n");
	    //inverse
	    for (int mat_row=0;mat_row<sizeC;mat_row++){
	      for (int mat_col=0;mat_col<sizeC;mat_col++){
		fprintf(pFile,"%12.4E",fjacinv(mat_row,mat_col));}
	      fprintf(pFile,"\n");} 
	    //solution
	    for (int mat_col=0;mat_col<sizeC;mat_col++){
	      fprintf(pFile,"%12.4E",pvec[mat_col]);} 
	      }	
	    */
	    
	    //check for convergence, update
	    T errx=0.;//%Check root convergence.
	    for (int row=0;row<sizeC;row++){
	      errx +=fabs(xvec[row]);
	      coeffg[c][row]+= xvec[row];
	    }
	    //if (c==_options.printCellNumber){cout <<"errx "<<errx<<endl;}
	    if(errx <= tolx)
	      break;
	    
	  }
	  
	}
      }
  }
  void EquilibriumDistributionESBGK()
  {
    ComputeMacroparametersESBGK();
    const int ktrial=_options.NewtonsMethod_ktrial;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();


	//const VectorT5Array& coeff = dynamic_cast<const VectorT5Array&>(_macroFields.coeff[cells]);
	//initialize coeffg
	VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
	/*
	for(int c=0; c<nCells;c++){
	  coeffg[c][0]=coeff[c][0];
	  coeffg[c][1]=coeff[c][1];
	  coeffg[c][2]=coeff[c][2];
	  coeffg[c][3]=coeff[c][1];
	  coeffg[c][4]=coeff[c][3];
	  coeffg[c][5]=coeff[c][1];
	  coeffg[c][6]=coeff[c][4];
	  coeffg[c][7]=0.0;
	  coeffg[c][8]=0.0;
	  coeffg[c][9]=0.0;	    
	}
	*/
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	NewtonsMethodESBGK(ktrial);
	
	for(int j=0;j< numFields;j++){
	  Field& fndEqES = *_dsfEqPtrES.dsf[j];
	  TArray& fEqES = dynamic_cast< TArray&>(fndEqES[cells]);
	  for(int c=0; c<nCells;c++){ 
	    T Cc1=(cx[j]-v[c][0]);
	    T Cc2=(cy[j]-v[c][1]);
	    T Cc3=(cz[j]-v[c][2]);
	    fEqES[c]=coeffg[c][0]*exp(-coeffg[c][1]*pow(Cc1,2)+coeffg[c][2]*Cc1
				      -coeffg[c][3]*pow(Cc2,2)+coeffg[c][4]*Cc2
				      -coeffg[c][5]*pow(Cc3,2)+coeffg[c][6]*Cc3
				    +coeffg[c][7]*cx[j]*cy[j]+coeffg[c][8]*cy[j]*cz[j]
				    +coeffg[c][9]*cz[j]*cx[j]);
	  }
	  
	}
      }
  }
  
  void setJacobianESBGK(SquareMatrix<T,10>& fjac, VectorT10& fvec, const VectorT10& xn,const VectorT3& v,const int c)
  {
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const TArray& wts = dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);	
    
    const TArray2D& malphaESBGK = dynamic_cast<const TArray2D&>(*_quadrature.malphaESBGKPtr);	
    const int numFields= _quadrature.getDirCount(); 
    VectorT10 mexp;
    
    for(int j=0;j< numFields;j++){   
      T Cc1=cx[j]-v[0];
      T Cc2=cy[j]-v[1];
      T Cc3=cz[j]-v[2];
      T Econst=xn[0]*exp(-xn[1]*pow(Cc1,2)+xn[2]*Cc1-xn[3]*pow(Cc2,2)+ xn[4]*Cc2
			 -xn[5]*pow(Cc3,2)+xn[6]*Cc3
			 +xn[7]*cx[j]*cy[j]+xn[8]*cy[j]*cz[j]+xn[9]*cz[j]*cx[j])*wts[j];     
      
      for (int row=0;row<10;row++){
	fvec[row]+= -Econst*malphaESBGK(j,row); //smm
      }
      
      mexp[0]=-Econst/xn[0];
      mexp[1]=Econst*pow(Cc1,2);
      mexp[2]=-Econst*Cc1;
      mexp[3]=Econst*pow(Cc2,2);
      mexp[4]=-Econst*Cc2; 
      mexp[5]=Econst*pow(Cc3,2);
      mexp[6]=-Econst*Cc3;
      
      mexp[7]=-Econst*cx[j]*cy[j];
      mexp[8]=-Econst*cy[j]*cz[j];
      mexp[9]=-Econst*cz[j]*cx[j];
   
      for (int row=0;row<10;row++){
	for (int col=0;col<10;col++){
	  fjac(row,col)+=malphaESBGK(j,row)*mexp[col];  //new
	}
      }
      
      
      
      
    }
    
        
  }

  void initializeMaxwellian()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();

	const double pi=_options.pi;
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 

	for(int j=0;j< numFields;j++){
	  Field& fnd = *_dsfPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=density[c]/pow((pi*temperature[c]),1.5)*
	      exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
		    pow((cz[j]-v[c][2]),2.0))/temperature[c]);
	    
	   
	  }
	  
	  if (_options.transient)
	    //updateTime();
	    
	    {
	      Field& fnd1 = *_dsfPtr1.dsf[j];
	      TArray& f1 = dynamic_cast< TArray&>(fnd1[cells]);
	      for (int c=0;c<nCells;c++)
		f1[c] = f[c];
              //cout << "discretization order " << _options.timeDiscretizationOrder << endl ;
	      if (_options.timeDiscretizationOrder > 1)
		{
		  Field& fnd2 = *_dsfPtr2.dsf[j];
		  TArray& f2 = dynamic_cast< TArray&>(fnd2[cells]);
		  for (int c=0;c<nCells;c++)
		    f2[c] = f[c];
		}
	    } 
	    
	  
	}
      }
  }

  void initializeFineMaxwellian()
  {
    const int numMeshes = _meshes.size();
    const T zero(0.);
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const int nCells = cells.getCount();

        const double pi=_options.pi;
        const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
        const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
        const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);

        const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
        const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
        const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
        const int numFields= _quadrature.getDirCount(); 

        for(int j=0;j< numFields;j++){
          Field& fnd0 = *_dsfPtr0.dsf[j];
	  Field& fndRes = *_dsfPtrRes.dsf[j];
          TArray& f0 = dynamic_cast< TArray&>(fnd0[cells]);
	  TArray& fRes = dynamic_cast< TArray&>(fndRes[cells]);
          for(int c=0; c<nCells;c++){
            f0[c]=density[c]/pow((pi*temperature[c]),1.5)*
              exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
                    pow((cz[j]-v[c][2]),2.0))/temperature[c]);
            fRes[c]=zero;
           
          }
        }
    }
  }

  void initializeCoarseMaxwellian()
  {
    const int numMeshes = _meshes.size();
    const T zero(0.);
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        const int nCells = cells.getCount();
        const int numFields= _quadrature.getDirCount();

        for(int j=0;j< numFields;j++){
          Field& fnd0 = *_dsfPtr0.dsf[j];
	  Field& fndFAS = *_dsfPtrFAS.dsf[j];
          Field& fndRes = *_dsfPtrRes.dsf[j];
          TArray& f0 = dynamic_cast< TArray&>(fnd0[cells]);
	  TArray& fFAS = dynamic_cast< TArray&>(fndFAS[cells]);
          TArray& fRes = dynamic_cast< TArray&>(fndRes[cells]);
          for(int c=0; c<nCells;c++){
            f0[c]=zero;
	    fFAS[c]=zero;
            fRes[c]=zero;

          }
        }
      }
  }

  void weightedMaxwellian(double weight1,double uvel1,double vvel1,double wvel1,double uvel2,double vvel2,double wvel2,double temp1,double temp2)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	//double pi(acos(-1.0));
	const double pi=_options.pi;
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	for(int j=0;j< numFields;j++){
	  Field& fnd = *_dsfPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=weight1*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-uvel1),2.0)+pow((cy[j]-vvel1),2.0)+pow((cz[j]-wvel1),2.0))/temp1)
	      +(1-weight1)*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-uvel2),2.0)+pow((cy[j]-vvel2),2.0)+pow((cz[j]-wvel2),2.0))/temp2);
	  }
	  if (_options.transient)
	    {
	      Field& fnd1 = *_dsfPtr1.dsf[j];
	      TArray& f1 = dynamic_cast< TArray&>(fnd1[cells]);
	      for (int c=0;c<nCells;c++)
		f1[c] = f[c];
              //cout << "discretization order " << _options.timeDiscretizationOrder << endl ;
	      if (_options.timeDiscretizationOrder > 1)
		{
		  Field& fnd2 = *_dsfPtr2.dsf[j];
		  TArray& f2 = dynamic_cast< TArray&>(fnd2[cells]);
		  for (int c=0;c<nCells;c++)
		    f2[c] = f[c];
		}
	    }
	}
      }
  }
  void weightedMaxwellian(double weight1,double uvel1,double uvel2,double temp1,double temp2)
  {
    const double vvel1=0.0;
    const double wvel1=0.0;
    const double vvel2=0.0;
    const double wvel2=0.0;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	//double pi(acos(-1.0));
	const double pi=_options.pi;
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	for(int j=0;j< numFields;j++){
	  Field& fnd = *_dsfPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=weight1*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-uvel1),2.0)+pow((cy[j]-vvel1),2.0)+pow((cz[j]-wvel1),2.0))/temp1)
	      +(1-weight1)*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-uvel2),2.0)+pow((cy[j]-vvel2),2.0)+pow((cz[j]-wvel2),2.0))/temp2);
	  }
	  if (_options.transient)
	    {
	      Field& fnd1 = *_dsfPtr1.dsf[j];
	      TArray& f1 = dynamic_cast< TArray&>(fnd1[cells]);
	      for (int c=0;c<nCells;c++)
		f1[c] = f[c];
              //cout << "discretization order " << _options.timeDiscretizationOrder << endl ;
	      if (_options.timeDiscretizationOrder > 1)
		{
		  Field& fnd2 = *_dsfPtr2.dsf[j];
		  TArray& f2 = dynamic_cast< TArray&>(fnd2[cells]);
		  for (int c=0;c<nCells;c++)
		    f2[c] = f[c];
		}
	    }
	}
      }
  }
  COMETBCMap& getBCMap() {return _bcMap;}
  COMETVCMap& getVCMap() {return _vcMap;}
  
  COMETModelOptions<T>&   getOptions() {return _options;}
  const map<int, vector<int> >&  getFaceReflectionArrayMap() const { return _faceReflectionArrayMap;}
  
 
  
  // const vector<int>& vecReflection = _faceReflectionArrayMap[faceID]
map<string,shared_ptr<ArrayBase> >&
  getPersistenceData()
  {
    _persistenceData.clear();
    
    Array<int>* niterArray = new Array<int>(1);
    (*niterArray)[0] = _niters;
    _persistenceData["niters"]=shared_ptr<ArrayBase>(niterArray);
    
    if (_initialKmodelNorm)
    {
      // _persistenceData["initialKmodelNorm"] =_initialKmodelNorm->getArrayPtr(_macroFields.pressure);
       const Field& dsfField = *_dsfPtr.dsf[0];
      _persistenceData["initialKmodelNorm"] =_initialKmodelNorm->getArrayPtr(dsfField);

    }
    else
    {
         Array<T>* xArray = new Array<T>(1);
        xArray->zero();
        _persistenceData["initialKmodelNorm"]=shared_ptr<ArrayBase>(xArray);
        
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

    if (_persistenceData.find("initialKmodelNorm") != _persistenceData.end())
    {
        shared_ptr<ArrayBase>  r = _persistenceData["initialKmodelNorm"];
        _initialKmodelNorm = MFRPtr(new MultiFieldReduction());
	Field& dsfField = *_dsfPtr.dsf[0];
        _initialKmodelNorm->addArray(dsfField,r);
	//_initialKmodelNorm->addArray(_dsfPtr,r);
    }
  }
  
  void SetBoundaryConditions()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	
	COMETVC<T> *vc(new COMETVC<T>());
	vc->vcType = "flow";
	_vcMap[mesh.getID()] = vc;
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if (_bcMap.find(fg.id) == _bcMap.end())
	      {
		COMETBC<T> *bc(new COMETBC<T>());
		
		_bcMap[fg.id] = bc;
		if((fg.groupType == "wall"))
		  {
		    bc->bcType = "WallBC";
		  }
		else if((fg.groupType == "realwall"))
		  {
		    bc->bcType = "RealWallBC";
		  }
		else if (fg.groupType == "velocity-inlet")
		  {
		    bc->bcType = "VelocityInletBC";
		  }
		else if (fg.groupType == "pressure-inlet")
		  {
		    bc->bcType = "PressureInletBC";
		  }
		else if (fg.groupType == "pressure-outlet")
		  {
		    bc->bcType = "PressureOutletBC";
		  }
		else if ((fg.groupType == "symmetry"))
		  {
		    bc->bcType = "SymmetryBC";
		    }
	
		else if((fg.groupType =="zero-gradient "))
		  {
		      bc->bcType = "ZeroGradBC";
		  }
		else
		  throw CException("COMETModel: unknown face group type "
				     + fg.groupType);
	      }
	  }
	/*
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if (_bcMap.find(fg.id) == _bcMap.end())
	      { COMETBC<T> *bc(new COMETBC<T>());
		
		_bcMap[fg.id] = bc;
		
		if ((fg.groupType == "NSinterface"))
		  {
		    bc->bcType = "NSInterfaceBC";
		  }
	      }
	  }
	*/
      }
  }
  
  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh = *_meshes[n];

	const StorageSite& cells = mesh.getCells();
	const int numFields= _quadrature.getDirCount(); 
	//callBoundaryConditions();  //new
	for (int direction = 0; direction < numFields; direction++)
	  {
	    Field& fnd = *_dsfPtr.dsf[direction];
	    Field& fndN1 = *_dsfPtr1.dsf[direction];
	    TArray& f = dynamic_cast<TArray&>(fnd[cells]);
	    TArray& fN1 = dynamic_cast<TArray&>(fndN1[cells]);
	    if (_options.timeDiscretizationOrder > 1)
	      {
		Field& fndN2 = *_dsfPtr2.dsf[direction];
		TArray& fN2 = dynamic_cast<TArray&>(fndN2[cells]);
		fN2 = fN1; 
	      }
	    fN1 = f;
	  }
#ifdef FVM_PARALLEL
	if ( MPI::COMM_WORLD.Get_rank() == 0 )
	  {cout << "updated time" <<endl;}
	
#endif
#ifndef FVM_PARALLEL 
	cout << "updated time" <<endl;
#endif
	//ComputeMacroparameters();	//update macroparameters
        //ComputeCollisionfrequency();
	//if (_options.fgamma==0){initializeMaxwellianEq();}
	//else{ EquilibriumDistributionBGK();}	
	//if (_options.fgamma==2){EquilibriumDistributionESBGK();}
	
      }
  }
	  

  void callCOMETBoundaryConditions()    
  { 
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            //const int nFaces = faces.getCount();
            const COMETBC<T>& bc = *_bcMap[fg.id];
            COMETBoundaryConditions<T,T,T> cbc(faces, mesh,_geomFields,_quadrature,_macroFields,_dsfPtr);
            FloatValEvaluator<T> bTemperature(bc.getVal("specifiedTemperature"),faces);
            FloatValEvaluator<T> bPressure(bc.getVal("specifiedPressure"),faces);
            if(bc.bcType=="PressureInletBC")
	      {
                cbc.applyPressureInletBC(bTemperature,bPressure);
	      } 
            else if(bc.bcType=="PressureOutletBC")
	      {
                cbc.applyPressureOutletBC(bTemperature,bPressure);
	      }
	  }
        foreach(const FaceGroupPtr igPtr, mesh.getInterfaceGroups())
	  {
            
            const FaceGroup& ig = *igPtr;
            const StorageSite& faces = ig.site;
            //const int nFaces = faces.getCount();
            
            
            COMETBoundaryConditions<T,T,T> cbc(faces, mesh,_geomFields,_quadrature,_macroFields,_dsfPtr); 
            if(ig.groupType=="NSinterface")
	      {
                cbc.applyNSInterfaceBC();//bTemperature,bPressure,bVelocity,bStress);
	      } 
	  }
      }//end of loop through meshes
  }
  
  void OutputDsfBLOCK(const char* filename)
  {
    FILE * pFile;
    pFile = fopen(filename,"w"); 
    int N1=_quadrature.getNVCount();
    int N2=_quadrature.getNthetaCount();
    int N3=_quadrature.getNphiCount();
    fprintf(pFile,"%s \n", "VARIABLES= cx, cy, cz, f,fEq,fES");
    fprintf(pFile, "%s %i %s %i %s %i \n","ZONE I=", N3,",J=",N2,",K=",N1);
    fprintf(pFile, "%s \n","F=BLOCK, VARLOCATION=(NODAL,NODAL,NODAL,NODAL,NODAL,NODAL)");
    const int numMeshes = _meshes.size();
    const int cellno=_options.printCellNumber;
    for (int n=0; n<numMeshes; n++){
      const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
      const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
      const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
      const int numFields= _quadrature.getDirCount(); 	
      for(int j=0;j< numFields;j++){fprintf(pFile,"%E \n",cx[j]);}
      for(int j=0;j< numFields;j++){fprintf(pFile,"%E \n",cy[j]);}
      for(int j=0;j< numFields;j++){fprintf(pFile,"%E \n",cz[j]);}
      
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      for(int j=0;j< numFields;j++){
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	fprintf(pFile,"%E\n",f[cellno]);
      } 
      for(int j=0;j< numFields;j++){
	Field& fEqnd = *_dsfEqPtr.dsf[j];
	TArray& fEq = dynamic_cast< TArray&>(fEqnd[cells]);
	fprintf(pFile,"%E\n",fEq[cellno]);
      }
      if(_options.fgamma==2){
	for(int j=0;j< numFields;j++){
	 
	    Field& fndEqES = *_dsfEqPtrES.dsf[j];
	  TArray& fEqES = dynamic_cast< TArray&>(fndEqES[cells]);
	  fprintf(pFile,"%E\n",fEqES[cellno]);
	}}
      else{
	for(int j=0;j< numFields;j++){

	    Field& fndEq = *_dsfEqPtr.dsf[j];
	  TArray& fEq = dynamic_cast< TArray&>(fndEq[cells]);
	  fprintf(pFile,"%E\n",fEq[cellno]);
	}}
      
    }
    fclose(pFile);
  }

  void MakeCoarseMesh1(const MeshList& inMeshes, GeomFields& inGeomFields,
                     MeshList& outMeshes)
  {

    int smallestMesh=-1;
    const int numMeshes=inMeshes.size();
    
    for(int n=0;n<numMeshes;n++)
      {
        const Mesh& mesh=*inMeshes[n];
        const int dim=mesh.getDimension();
        Mesh* newMeshPtr=new Mesh(dim);

        outMeshes.push_back(newMeshPtr);
	
	const StorageSite& inCells=mesh.getCells();
        StorageSite& outCells=newMeshPtr->getCells();
        StorageSite& outFaces=newMeshPtr->getFaces();
        const StorageSite& inFaces=mesh.getFaces();
        const int inCellCount=inCells.getSelfCount();
        const int inCellTotal=inCells.getCount();
        const int inFaceCount=inFaces.getCount();
        const int inGhost=inCellTotal-inCellCount;
        int coarseCount=0;
        _siteMap[&inCells]=&outCells;
        Field& FineToCoarseField=inGeomFields.fineToCoarse;
        IntArray& FineToCoarse=dynamic_cast<IntArray&>(FineToCoarseField[inCells]);

        const CRConnectivity& inCellinFaces=mesh.getCellFaces();
        const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
        Field& areaMagField=inGeomFields.areaMag;
        Field& areaField=inGeomFields.area;
        const TArray& areaMagArray=dynamic_cast<const TArray&>(areaMagField[inFaces]);
	const VectorT3Array& areaArray=dynamic_cast<const VectorT3Array&>(areaField[inFaces]);
        const BCfaceArray& inBCfArray=*(_BFaces[n]);

	const IntArray& ibType = dynamic_cast<const IntArray&>(inGeomFields.ibType[inCells]);

	//first sweep to make initial pairing
        int pairWith;
        Array<bool> marker(inFaceCount);
        Array<bool> marked(inCellCount);
	const T zero(0.);
	for(int c=0;c<inCellCount;c++)
        {
	    if((FineToCoarse[c]<0)&&(ibType[c] == Mesh::IBTYPE_FLUID)) //dont bother if im already paired
            {
                //loop through all neighbors to find pairing
                const int neibCount=inCellinFaces.getCount(c);
                pairWith=-1;
                T maxArea=0.;
                int c2;
                for(int i=0; i<inFaceCount; i++)
		{
                    marker[i] = false;
		}
                for(int i=0; i<inCellCount; i++)
		{
                    marked[i] = false;
		}
                for(int neib=0;neib<neibCount;neib++)
                {
                    const int f=inCellinFaces(c,neib);
                    
                    if(inBCfArray[f]==0)  //not a boundary face
                    {
                        if(c==inFaceinCells(f,1))
                          c2=inFaceinCells(f,0);
                        else
                          c2=inFaceinCells(f,1);

                        VectorT3 tempArea;
			tempArea[0]=zero;
			tempArea[1]=zero;
			tempArea[2]=zero;
                        if((FineToCoarse[c2]==-1)&&(!marked[c2]))
			{
                            marker[f]=true;
                            marked[c2]=true;
                            for(int face=0;face<inFaceCount;face++)
			    {
                                if((c==inFaceinCells(face,0))&&(c2==inFaceinCells(face,1)))
                                  tempArea+=areaArray[face];
                                else if((c2==inFaceinCells(face,0))&&(c==inFaceinCells(face,1)))
                                  tempArea-=areaArray[face];
			    }
			}

                        if((FineToCoarse[c2]==-1)&&(marker[f]))
                          if(mag(tempArea)>maxArea)
			  {
                              pairWith=c2;
                              maxArea=mag(tempArea);
			  }

			/*
                        if(FineToCoarse[c2]==-1)
                          if(areaMagArray[f]>maxArea)
			  {
			      pairWith=c2;
			      maxArea=areaMagArray[f];
			  }
			*/
		    }
		}
		
                if(pairWith!=-1)
                {
                    FineToCoarse[c]=coarseCount;
                    FineToCoarse[pairWith]=coarseCount;
                    coarseCount++;
		}
	    }
	}
	
        //second sweep to group stragglers, or group with self
        for(int c=0;c<inCellCount;c++)
        {
	    if((FineToCoarse[c]==-1)&&(ibType[c] == Mesh::IBTYPE_FLUID))
            {
                const int neibCount=inCellinFaces.getCount(c);
                T maxArea=0.;
                int c2,c2perm;
                pairWith=-2;
                for(int i=0; i<inFaceCount; i++)
		{
                    marker[i] = false;
		}
                for(int i=0; i<inCellCount; i++)
		{
                    marked[i] = false;
		}

                for(int neib=0;neib<neibCount;neib++)
                {
                    const int f=inCellinFaces(c,neib);
                    
                    if(inBCfArray[f]==0)  //not a boundary face
                    {
                        if(c==inFaceinCells(f,1))
                          c2=inFaceinCells(f,0);
                        else
                          c2=inFaceinCells(f,1);

			VectorT3 tempArea;
                        tempArea[0]=zero;
                        tempArea[1]=zero;
                        tempArea[2]=zero;
                        if(!marked[c2])
			{
                            marker[f]=true;
                            marked[c2]=true;
                            for(int face=0;face<inFaceCount;face++)
			    {
                                if((c==inFaceinCells(face,0))&&(c2==inFaceinCells(face,1)))
                                  tempArea+=areaArray[face];
                                else if((c2==inFaceinCells(face,0))&&(c==inFaceinCells(face,1)))
                                  tempArea-=areaArray[face];
			    }
			}
			
                        if(marker[f])
                          if(mag(tempArea)>maxArea)
			  {
                              pairWith=FineToCoarse[c2]; //coarse level cell
                              c2perm=c2;                 //fine level cell
                              maxArea=mag(tempArea);
			  }

			/*
                        if(areaMagArray[f]>maxArea)
                        {
                            pairWith=FineToCoarse[c2]; //coarse level cell
                            c2perm=c2;                 //fine level cell
			    maxArea=areaMagArray[f];
			}
			*/
		    }
		}

                if(pairWith==-2)
                {
                    FineToCoarse[c]=coarseCount;
                    coarseCount++;
		}
                else
                {
                    if(FineToCoarse[c2perm]==-1)
                    {
                        FineToCoarse[c]=coarseCount;
                        FineToCoarse[c2perm]=coarseCount;
                        coarseCount++;
		    }
                    else
                      FineToCoarse[c]=pairWith;
		}
	    }     
	}

        for(int c=0;c<inCellCount;c++)
	{
            if(ibType[c] != Mesh::IBTYPE_FLUID)
	    {
                FineToCoarse[c]=coarseCount;
                coarseCount++;
	    }
	}

        int coarseGhost=coarseCount;
	_coarseSizes[&mesh]=coarseCount;

	/*	
	if(MPI::COMM_WORLD.Get_rank()==1)
	  for(int c=0;c<inCells.getCount();c++)
	    cout<<" in makecoarsemesh1 before boundary, rank,level,cell,index = "<<MPI::COMM_WORLD.Get_rank()<<" "<<_level<<" "<<c<<" "<<FineToCoarse[c]<<endl;
	*/	

        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            const int nFaces = faces.getCount();

            const CRConnectivity& faceCells = mesh.getFaceCells(faces);                  
            
	    for(int f=0; f< nFaces; f++)
	    {
		const int c1= faceCells(f,1);// boundary cell
		FineToCoarse[c1]=coarseGhost;
		coarseGhost++;
	    }
	}
	
	/*
	if(MPI::COMM_WORLD.Get_rank()==1)
	  for(int c=0;c<inCells.getCount();c++)
	    cout<<" in makecoarsemesh1, rank, level,cell,index = "<<MPI::COMM_WORLD.Get_rank()<<" "<<_level<<" "<<c<<" "<<FineToCoarse[c]<<endl;
	*/
      }
  }

  int MakeCoarseMesh2(const MeshList& inMeshes, GeomFields& inGeomFields,
		      GeomFields& coarseGeomFields,MeshList& outMeshes)
  {

    int smallestMesh=-1;
    const int numMeshes=inMeshes.size();
    for(int n=0;n<numMeshes;n++)
      {
        const Mesh& mesh=*inMeshes[n];
        const int dim=mesh.getDimension();
        //Mesh* newMeshPtr=new Mesh(dim);

        //outMeshes.push_back(newMeshPtr);
	Mesh& newMeshPtr=*outMeshes[n];

        const StorageSite& inCells=mesh.getCells();
        StorageSite& outCells=newMeshPtr.getCells();
        StorageSite& outFaces=newMeshPtr.getFaces();
	StorageSite& outIBFaces=newMeshPtr.getIBFaces();
	const IntArray& ibType = dynamic_cast<const IntArray&>(inGeomFields.ibType[inCells]);
        const StorageSite& inFaces=mesh.getFaces();
	const StorageSite& inIBFaces=mesh.getIBFaces();
        const int inCellCount=inCells.getSelfCount();
        const int inCellTotal=inCells.getCount();
        const int inFaceCount=inFaces.getCount();
        const int inGhost=inCellTotal-inCellCount;
	int outGhost = 0;
        int coarseCount=0;
        Field& FineToCoarseField=inGeomFields.fineToCoarse;
	const IntArray& FineToCoarse=dynamic_cast<const IntArray&>(FineToCoarseField[inCells]);
	const CRConnectivity& inCellinFaces=mesh.getCellFaces();
        const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
        Field& areaMagField=inGeomFields.areaMag;
        const TArray& areaMagArray=dynamic_cast<const TArray&>(areaMagField[inFaces]);
        const BCfaceArray& inBCfArray=*(_BFaces[n]);

	/*** creating coarse level cells ***/
	coarseCount = _coarseSizes[&mesh];
	int interfaceCellsLevel0 = _coarseGhostSizes[&mesh];	

	int boundaryCell=0;
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
	    
	    const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	    
	    for(int f=0; f< nFaces; f++)
	      {
		const int c1= faceCells(f,1);// boundary cell
		if(boundaryCell<FineToCoarse[c1])
		  boundaryCell=FineToCoarse[c1];
	      }
	}
	boundaryCell++;
	if(boundaryCell!=1)
	  boundaryCell-=coarseCount;
	else
	  boundaryCell=0;

	for(int c=0; c< inCells.getCountLevel1(); c++)
	{
            if(outGhost<FineToCoarse[c])
              outGhost=FineToCoarse[c];
	}
        outGhost++;
        outGhost-=coarseCount;

	int interfaceCells = outGhost - boundaryCell;

	/*	
	if(MPI::COMM_WORLD.Get_rank()==1)
	  cout<<"rank,level,inCellinternal, incelltotal,outCellInternal, outCellExternal, outInterface ="<<MPI::COMM_WORLD.Get_rank()<<" "<<_level<<" "<<inCellCount<<" and  "<<inCellTotal<<" and "<<coarseCount<<" and "<<outGhost<<" "<<interfaceCells<<endl;
	*/
	
	/*
	if(MPI::COMM_WORLD.Get_rank()==1)
	  for(int c=0;c<inCellTotal;c++)
	    cout<<" in makecoarsemesh2, rank,level, cell,index = "<<MPI::COMM_WORLD.Get_rank()<<" "<<_level<<" "<<c<<" and "<<FineToCoarse[c]<<endl;
	*/

        outCells.setCount(coarseCount,boundaryCell+interfaceCellsLevel0);
        outCells.setCountLevel1(coarseCount+outGhost);        

	/*** created coarse level cells ***/

	//make the coarse cell to fine cell connectivity.
	CRPtr CoarseToFineCells=CRPtr(new CRConnectivity(outCells,inCells));
        CoarseToFineCells->initCount();
	
        for(int c=0;c<inCells.getCountLevel1();c++)
          CoarseToFineCells->addCount(FineToCoarse[c],1);

        CoarseToFineCells->finishCount();

        for(int c=0;c<inCells.getCountLevel1();c++)
          CoarseToFineCells->add(FineToCoarse[c],c);

        CoarseToFineCells->finishAdd();

        /*** connectivity between itself (cells) and its finer mesh cells ***/ 
        newMeshPtr.setConnectivity(outCells,inCells,CoarseToFineCells);

        CRPtr FineFacesCoarseCells=CRPtr(new CRConnectivity(inFaces,outCells));
        FineFacesCoarseCells->initCount();

        //count surviving faces
        int survivingFaces=0;
        int coarse0, coarse1;
        for(int f=0;f<inFaceCount;f++)
        {
            coarse0=FineToCoarse[inFaceinCells(f,0)];
            coarse1=FineToCoarse[inFaceinCells(f,1)];
            if(coarse0!=coarse1)
            {
                survivingFaces++;
                FineFacesCoarseCells->addCount(f,2);
	    }
	}
	
        FineFacesCoarseCells->finishCount();

        //make non-zero's
        int fc0,fc1,cc0,cc1;
        for(int f=0;f<inFaceCount;f++)
          {
            fc0=inFaceinCells(f,0);
            fc1=inFaceinCells(f,1);
            cc0=FineToCoarse[fc0];
            cc1=FineToCoarse[fc1];
            if(cc0!=cc1)
              {
                FineFacesCoarseCells->add(f,cc0);
                FineFacesCoarseCells->add(f,cc1);
              }
          }
	
        FineFacesCoarseCells->finishAdd();

        CRPtr CoarseCellsFineFaces=FineFacesCoarseCells->getTranspose();
        CRPtr CellCellCoarse=CoarseCellsFineFaces->multiply(*FineFacesCoarseCells,true);

	/*** coarse level cellcell connectivity created ***/ 

	/*
        int counter=0;
        BArray counted(outCells.getCount());
        counted=false;
        for(int c=0;c<outCells.getCount();c++)
          {
            counted[c]=true;
            const int neibs=CellCellCoarse->getCount(c);
            for(int n=0;n<neibs;n++)
              {
                const int c1=(*CellCellCoarse)(c,n);
                if(!counted[c1])
                  counter++;
              }
          }

	//outFaces.setCount(counter);
	*/

        int countFaces=0;
        int cCell0, cCell1;
        for(int f=0;f<inFaceCount;f++)
        {
            cCell0=FineToCoarse[inFaceinCells(f,0)];
            cCell1=FineToCoarse[inFaceinCells(f,1)];
            if(cCell0!=cCell1)
            {
                countFaces++;
	    }
	}

	outFaces.setCount(countFaces);

	int inFaceGhost = 0;
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  inFaceGhost+=(*fgPtr).site.getCount();
        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  inFaceGhost+=(*fgPtr).site.getCount();
	const int inInteriorFaces = inFaceCount - inFaceGhost;

	const int del = inFaceCount - outFaces.getCount();

	//const int interiorCount=outFaces.getCount()-inGhost;
	const int interiorCount=inInteriorFaces-del;
	//if(MPI::COMM_WORLD.Get_rank()==1)
	//cout<<"level,outfaces, outghost, totalfaces = "<<_level<<" "<<interiorCount<<" "<<inGhost<<" "<<countFaces<<endl;
	const StorageSite& interiorFaces=newMeshPtr.createInteriorFaceGroup(interiorCount);
	
        int inOffset=interiorCount;
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	{
            const FaceGroup& fg=*fgPtr;
            const int size=fg.site.getCount();
            newMeshPtr.createBoundaryFaceGroup(size,inOffset,fg.id,fg.groupType);
            inOffset+=size;
	}

        foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	{
            const FaceGroup& fg=*fgPtr;
            const int size=fg.site.getCount();
            newMeshPtr.createInterfaceGroup(size,inOffset,fg.id);
            inOffset+=size;
	}

        CRPtr CoarseFaceCoarseCell=CRPtr(new CRConnectivity(outFaces,outCells));
        CoarseFaceCoarseCell->initCount();

        survivingFaces=0;
        for(int f=0;f<inFaceCount;f++)
        {
            coarse0=FineToCoarse[inFaceinCells(f,0)];
            coarse1=FineToCoarse[inFaceinCells(f,1)];
            if(coarse0!=coarse1)
            {
                CoarseFaceCoarseCell->addCount(survivingFaces,2);
		survivingFaces++;
	    }
	}
	
        CoarseFaceCoarseCell->finishCount();

        //make non-zero's
	survivingFaces=0;
        for(int f=0;f<inFaceCount;f++)
          {
            fc0=inFaceinCells(f,0);
            fc1=inFaceinCells(f,1);
            cc0=FineToCoarse[fc0];
            cc1=FineToCoarse[fc1];
            if(cc0!=cc1)
              {
                CoarseFaceCoarseCell->add(survivingFaces,cc0);
                CoarseFaceCoarseCell->add(survivingFaces,cc1);
		survivingFaces++;
              }
          }
	
        CoarseFaceCoarseCell->finishAdd();

	CRPtr CoarseCellCoarseFace=CoarseFaceCoarseCell->getTranspose();

        newMeshPtr.setConnectivity(outCells,outFaces,CoarseCellCoarseFace);
        newMeshPtr.setConnectivity(outFaces,outCells,CoarseFaceCoarseCell);

        Field& coarseIbTypeField=coarseGeomFields.ibType;
        shared_ptr<IntArray> ibTypePtr(new IntArray(outCells.getCountLevel1()));
        *ibTypePtr = Mesh::IBTYPE_FLUID;
        coarseIbTypeField.addArray(outCells,ibTypePtr);

        IntArray& coarseIBType = dynamic_cast<IntArray&>(coarseGeomFields.ibType[outCells]);

        for(int c=0;c<inCellCount;c++)
          if(ibType[c] != Mesh::IBTYPE_FLUID)
            coarseIBType[FineToCoarse[c]]=ibType[c];

        outIBFaces.setCount(inIBFaces.getCount());
        _siteMap[&inIBFaces]=&outIBFaces;

        shared_ptr<IntArray> ibFaceIndexPtr(new IntArray(outFaces.getCount()));
        *ibFaceIndexPtr = -1;
        coarseGeomFields.ibFaceIndex.addArray(outFaces,ibFaceIndexPtr);

        const IntArray& fineIBFaceIndex=dynamic_cast<const IntArray&>(inGeomFields.ibFaceIndex[inFaces]);
        IntArray& coarseIBFaceIndex=dynamic_cast<IntArray&>(coarseGeomFields.ibFaceIndex[outFaces]);

        CRPtr CoarseFacesFineFaces=CRPtr(new CRConnectivity(outFaces,inFaces));
        CoarseFacesFineFaces->initCount();

	survivingFaces=0;
	for(int f=0;f<inFaceCount;f++)
          {
            int fc0=inFaceinCells(f,0);
            int fc1=inFaceinCells(f,1);
            const int cc0=FineToCoarse[fc0];
            const int cc1=FineToCoarse[fc1];

            if(cc1!=cc0)
              {
		CoarseFacesFineFaces->addCount(survivingFaces,1);
		survivingFaces++;
              }
          }
	//cout<<"hello 3 rank "<<MPI::COMM_WORLD.Get_rank()<<endl;
        CoarseFacesFineFaces->finishCount();
	
	survivingFaces=0;
	for(int f=0;f<inFaceCount;f++)
        {
            int fc0=inFaceinCells(f,0);
            int fc1=inFaceinCells(f,1);
            const int cc0=FineToCoarse[fc0];
            const int cc1=FineToCoarse[fc1];
            if(cc1!=cc0)
            {
		CoarseFacesFineFaces->add(survivingFaces,f);
		coarseIBFaceIndex[survivingFaces]=fineIBFaceIndex[f];
		survivingFaces++;
	    }
	}
	//cout<<"hello 4 rank "<<MPI::COMM_WORLD.Get_rank()<<endl;
        CoarseFacesFineFaces->finishAdd();

	/*
	if(MPI::COMM_WORLD.Get_rank()==1)
	{
	    for(int f=0;f<outFaces.getCount();f++)
	    {
		cout<<"level,rank,face,neighbor = "<<_level<<" "<<MPI::COMM_WORLD.Get_rank()<<" "<<f<<" "<<(*CoarseFaceCoarseCell)(f,0)<<" "<<(*CoarseFaceCoarseCell)(f,1)<<endl;
	    }
	}
	*/	

	//cout<<"hello 5 rank "<<MPI::COMM_WORLD.Get_rank()<<endl;
	
        //now make the geom fields
        const int outCellsCount=outCells.getSelfCount();
        TArrptr outCellVolumePtr=TArrptr(new TArray(outCellsCount));
        TArray& outCV=*outCellVolumePtr;
        outCV=0.;
	
	Field& VolumeField=inGeomFields.volume;
	Field& coarseVolumeField=coarseGeomFields.volume;
        const TArray& inCV=dynamic_cast<const TArray&>(VolumeField[inCells]);

        for(int c=0;c<outCellsCount;c++)
          {
            const int fineCount=CoarseToFineCells->getCount(c);
            for(int i=0;i<fineCount;i++)
              {
                int fc=(*CoarseToFineCells)(c,i);
                outCV[c]+=inCV[fc];
              }
          }

        coarseVolumeField.addArray(outCells,outCellVolumePtr);

	//cout<<"hello 6"<<endl;
	
	const int outFacesCount=outFaces.getCount();
        VT3Ptr outFaceAreaPtr=VT3Ptr(new VectorT3Array(outFacesCount));
        VectorT3Array& outFA=*outFaceAreaPtr;
        TArrptr outFaceAreaMagPtr=TArrptr(new TArray(outFacesCount));
        TArray& outFAMag=*outFaceAreaMagPtr;

        Field& FaceAreaField=inGeomFields.area;
	Field& coarseFaceAreaField=coarseGeomFields.area;
	Field& coarseAreaMagField=coarseGeomFields.areaMag;
        const VectorT3Array& inFA=
          dynamic_cast<const VectorT3Array&>(FaceAreaField[inFaces]);

        VectorT3 myZero;
        myZero[0]=0.;
        myZero[1]=0.;
        myZero[2]=0.;

        outFA=myZero;
        outFAMag=0.;
        for(int f=0;f<outFacesCount;f++)
          {
            const int fineCount=CoarseFacesFineFaces->getCount(f);
            const int cCell0=(*CoarseFaceCoarseCell)(f,0);
            for(int i=0;i<fineCount;i++)
              {
                const int fFace=(*CoarseFacesFineFaces)(f,i);
                const int fCell0=inFaceinCells(fFace,0);
                const int CCell0=FineToCoarse[fCell0];

                //must make sure the area vector is pointing
                //from c0 to c1
                if(CCell0==cCell0)
                  outFA[f]+=inFA[fFace];
                else
                  outFA[f]-=inFA[fFace];
                outFAMag[f]+=areaMagArray[fFace];
              }
          }
	
	coarseFaceAreaField.addArray(outFaces,outFaceAreaPtr);
        coarseAreaMagField.addArray(outFaces,outFaceAreaMagPtr);
	
	//cout<<"hello 7"<<endl;
	/*
        Field& ibTypeField=inGeomFields.ibType;
	Field& coarseIbTypeField=coarseGeomFields.ibType;
        shared_ptr<IntArray> ibTypePtr(new IntArray(outCells.getSelfCount()));
        *ibTypePtr = Mesh::IBTYPE_FLUID;
        coarseIbTypeField.addArray(outCells,ibTypePtr);
	*/
        if(smallestMesh<0)
          smallestMesh=outCells.getSelfCount();
        else
          {
            if(outCells.getSelfCount()<smallestMesh)
              smallestMesh=outCells.getSelfCount();
          }
	//cout<<"hello 8"<<endl;
	/*
        //This is for checking purposes only
        cout<<"Coarse Faces to Fine Faces"<<endl;
        for(int f=0;f<outFaces.getCount();f++)
          {
            const int neibs=CoarseFacesFineFaces->getCount(f);
            for(int n=0;n<neibs;n++)
              cout<<f<<" "<<(*CoarseFacesFineFaces)(f,n)<<endl;
            cout<<endl;
          }
        cout<<"Coarse Cells to Coarse Faces"<<endl;
        for(int c=0;c<outCells.getCount();c++)
          {
            const int neibs=CoarseCellCoarseFace->getCount(c);
            for(int n=0;n<neibs;n++)
              cout<<c<<" "<<(*CoarseCellCoarseFace)(c,n)<<endl;
            cout<<endl;
          }
	*/
      }
    //cout<<"hello 6"<<endl;
    return smallestMesh;
  }

  
  void syncGhostCoarsening(const MeshList& inMeshes, GeomFields& inGeomFields,
			   MeshList& outMeshes)
  {
  //const int xLen = coarseIndexField.getLength();

  //#pragma omp parallel for
    const int numMeshes=inMeshes.size();
    for(int n=0;n<numMeshes;n++)
    {
      const Mesh& mesh=*inMeshes[n];
      const int dim=mesh.getDimension();
      Mesh& newMeshPtr=*outMeshes[n];

      const StorageSite& inCells=mesh.getCells();
      const StorageSite& site=mesh.getCells();
      //const StorageSite& site=newMeshPtr.getCells();
      const int inCellCount=inCells.getSelfCount();
      const int inCellTotal=inCells.getCount();
      
      Field& FineToCoarseField=inGeomFields.fineToCoarse;
      IntArray& coarseIndex=dynamic_cast<IntArray&>(FineToCoarseField[inCells]);
      IntArray tempIndex(inCells.getCountLevel1());
      for(int c=0;c<inCells.getCountLevel1();c++)
        tempIndex[c]=coarseIndex[c];

      int coarseGhostSize=0;
      int tempGhostSize=0;
      //const int coarseSize = _coarseSizes.find(rowIndex)->second;
      //const int coarseSize = _coarseSizes[&mesh];
      int coarseSize = -1;
      coarseSize = _coarseSizes[&mesh];

      int boundaryCell=0;
      foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	  const FaceGroup& fg = *fgPtr;
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();

	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	  for(int f=0; f< nFaces; f++)
          {
	      const int c1= faceCells(f,1);// boundary cell
	      //if(MPI::COMM_WORLD.Get_rank()==1)
	      //cout<<"fgid, boundarycell, boundary coarse index = "<<fg.id<<" "<<c1<<" "<<coarseIndex[c1]<<endl;
	      if(boundaryCell<coarseIndex[c1])
		boundaryCell=coarseIndex[c1];
	  }
      }
      boundaryCell++;
      if(boundaryCell!=1)
	coarseSize = boundaryCell;
  
      //const StorageSite& site = *rowIndex.second;

      const StorageSite::GatherMap&  gatherMap  = site.getGatherMap();
      const StorageSite::GatherMap&  gatherMapLevel1  = site.getGatherMapLevel1();

      // collect all the toIndices arrays for each storage site from
      // both gatherMap and gatherMapLevel1
      
      typedef map<const StorageSite*, vector<const Array<int>* > > IndicesMap;
      IndicesMap toIndicesMap;
      IndicesMap tempIndicesMap;

      foreach(const StorageSite::GatherMap::value_type pos, gatherMap)
	{
          const StorageSite& oSite = *pos.first;
	  const Array<int>& tempIndices = *pos.second;
          const Array<int>& toIndices = *pos.second;

	  tempIndicesMap[&oSite].push_back(&tempIndices);
          toIndicesMap[&oSite].push_back(&toIndices);
	}
      
      foreach(const StorageSite::GatherMap::value_type pos, gatherMapLevel1)
	{
          const StorageSite& oSite = *pos.first;
          const Array<int>& toIndices = *pos.second;

          int found=0;
          foreach(const StorageSite::GatherMap::value_type posLevel0, gatherMap)
	  {
              const StorageSite& oSiteLevel0 = *posLevel0.first;
              if(oSite.getTag()==oSiteLevel0.getTag())
	      {
                  toIndicesMap[&oSiteLevel0].push_back(&toIndices);
                  found=1;
	      }
	  }
          
          if(found==0)
            toIndicesMap[&oSite].push_back(&toIndices);
	}

      foreach(IndicesMap::value_type pos, tempIndicesMap)
      {
          const StorageSite& oSite = *pos.first;
          const vector<const Array<int>* > tempIndicesArrays = pos.second;


          map<int,int> otherToMyMapping;
          UnorderedSet  gatherSet;

          foreach(const Array<int>* tempIndicesPtr, tempIndicesArrays)
	  {
              const Array<int>& tempIndices = *tempIndicesPtr;
              const int nGhostRows = tempIndices.getLength();
              for(int ng=0; ng<nGhostRows; ng++)
	      {
                  const int fineIndex = tempIndices[ng];
                  const int coarseOtherIndex = tempIndex[fineIndex];

                  if (coarseOtherIndex < 0)
                    continue;


                  if (otherToMyMapping.find(coarseOtherIndex) !=
                      otherToMyMapping.end())
		  {
                      
                      tempIndex[fineIndex] = otherToMyMapping[coarseOtherIndex];
		  }
                  else
		  {
                      tempIndex[fineIndex] = tempGhostSize+coarseSize;
                      otherToMyMapping[coarseOtherIndex] = tempIndex[fineIndex];
                      gatherSet.insert( tempIndex[fineIndex] );
                      tempGhostSize++;
		  }
	      }
	  }
      }
      
      foreach(IndicesMap::value_type pos, toIndicesMap)
        {
          const StorageSite& oSite = *pos.first;
          const vector<const Array<int>* > toIndicesArrays = pos.second;

          
          map<int,int> otherToMyMapping;
          UnorderedSet  gatherSet;

          foreach(const Array<int>* toIndicesPtr, toIndicesArrays)
	    {
              const Array<int>& toIndices = *toIndicesPtr;
              const int nGhostRows = toIndices.getLength();
              for(int ng=0; ng<nGhostRows; ng++)
		{
                  const int fineIndex = toIndices[ng];
                  const int coarseOtherIndex = coarseIndex[fineIndex];
                  
                  if (coarseOtherIndex < 0)
                    continue;
                  
                  if (otherToMyMapping.find(coarseOtherIndex) !=
                      otherToMyMapping.end())
		    {
                      coarseIndex[fineIndex] = otherToMyMapping[coarseOtherIndex];
		    }
                  else
		    {
                      coarseIndex[fineIndex] = coarseGhostSize+coarseSize;
		      otherToMyMapping[coarseOtherIndex] = coarseIndex[fineIndex];
                      gatherSet.insert( coarseIndex[fineIndex] );
                      coarseGhostSize++;
		    }
		  //if(MPI::COMM_WORLD.Get_rank()==1)
		  //cout<<"level,fineIndex, coarseIndex = "<<_level<<" "<<fineIndex<<" "<<coarseIndex[fineIndex]<<endl;
		}
	      //if(MPI::COMM_WORLD.Get_rank()==1)
	      //cout<<endl<<endl;
	    }

          const int coarseMappersSize = otherToMyMapping.size();

          shared_ptr<Array<int> > coarseToIndices(new Array<int>(coarseMappersSize));

          for(int n = 0; n < gatherSet.size(); n++)
	    {
              (*coarseToIndices)[n]   = gatherSet.getData().at(n);
	    }

          SSPair sskey(&site,&oSite);
          _coarseGatherMaps [sskey] = coarseToIndices;
        }

      const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      const StorageSite::ScatterMap& scatterMapLevel1 = site.getScatterMapLevel1();
      
      IndicesMap fromIndicesMap;

      foreach(const StorageSite::GatherMap::value_type pos, scatterMap)
	{
          const StorageSite& oSite = *pos.first;
          const Array<int>& fromIndices = *pos.second;

          fromIndicesMap[&oSite].push_back(&fromIndices);
	}
      
      foreach(const StorageSite::GatherMap::value_type pos, scatterMapLevel1)
	{
          const StorageSite& oSite = *pos.first;
          const Array<int>& fromIndices = *pos.second;

          int found=0;
          foreach(const StorageSite::ScatterMap::value_type posLevel0, scatterMap)
	    {
              const StorageSite& oSiteLevel0 = *posLevel0.first;
              if(oSite.getTag()==oSiteLevel0.getTag())
		{
                  fromIndicesMap[&oSiteLevel0].push_back(&fromIndices);
                  found=1;
		}
	    }
          
          if(found==0)
            fromIndicesMap[&oSite].push_back(&fromIndices);
	}
      
      foreach(IndicesMap::value_type pos, fromIndicesMap)
	{
          const StorageSite& oSite = *pos.first;
          const vector<const Array<int>* > fromIndicesArrays = pos.second;

          UnorderedSet  scatterSet;

          foreach(const Array<int>* fromIndicesPtr, fromIndicesArrays)
	    {
              const Array<int>& fromIndices = *fromIndicesPtr;
              const int nGhostRows = fromIndices.getLength();
              for(int ng=0; ng<nGhostRows; ng++)
		{
                  const int fineIndex = fromIndices[ng];
                  const int coarseOtherIndex = coarseIndex[fineIndex];
                  if (coarseOtherIndex >= 0)
                    scatterSet.insert( coarseOtherIndex );
		}

	    }

          const int coarseMappersSize = scatterSet.size();
          
          shared_ptr<Array<int> > coarseFromIndices(new Array<int>(coarseMappersSize));
          
          for(int n = 0; n < scatterSet.size(); n++ ) {
	    (*coarseFromIndices)[n] = scatterSet.getData().at(n);
          }
          
          SSPair sskey(&site,&oSite);
          _coarseScatterMaps[sskey] = coarseFromIndices;
          
	}
      //_coarseGhostSizes[rowIndex]=coarseGhostSize;
      _coarseGhostSizes[&mesh]=coarseGhostSize;
    }
  }

 void computeIBFaceDsf(const StorageSite& solidFaces,const int method,const int RelaxDistribution=0)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    typedef CRMatrixTranspose<T,VectorT3,VectorT3> IMatrixV3;
    if (method==1){
    const int numMeshes = _meshes.size();
    const int numFields= _quadrature.getDirCount(); 
    for (int direction = 0; direction < numFields; direction++)
      {
	Field& fnd = *_dsfPtr.dsf[direction];
	const TArray& pV =
	  dynamic_cast<const TArray&>(fnd[solidFaces]);
 #ifdef FVM_PARALLEL
      	MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE,pV.getData(),solidFaces.getCount() , MPI::DOUBLE, MPI::SUM); 
 #endif 

	for (int n=0; n<numMeshes; n++)
	  {	    
	    const Mesh& mesh = *_meshes[n];
	    const Mesh& fMesh = *_finestMeshes[n];
	    if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	      
	      const StorageSite& cells = mesh.getCells();
	      const StorageSite& fCells = fMesh.getCells();
	      const StorageSite& ibFaces = mesh.getIBFaces();
	      const StorageSite& fIbFaces = fMesh.getIBFaces();
	      
	      GeomFields::SSPair key1(&fIbFaces,&fCells);
	      const IMatrix& mIC =
		dynamic_cast<const IMatrix&>
		(*_finestGeomFields._interpolationMatrices[key1]);
	      
	      IMatrix mICV(mIC);
	   

           GeomFields::SSPair key2(&fIbFaces,&solidFaces);
           const IMatrix& mIP =
	     dynamic_cast<const IMatrix&>
	     (*_finestGeomFields._interpolationMatrices[key2]);

           IMatrix mIPV(mIP);
	   

           shared_ptr<TArray> ibV(new TArray(ibFaces.getCount()));
 
	   Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
	   const Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<const Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);
      
           const TArray& cV =
            dynamic_cast<const TArray&>(fnd[cells]);
	   TArray cFV(fCells.getCountLevel1());
	   for(int c=0;c<fCells.getCountLevel1();c++)
	     cFV[c]=cV[FinestToCoarse[c][_level]];

           ibV->zero();

           mICV.multiplyAndAdd(*ibV,cFV);
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

          fnd.addArray(ibFaces,ibV);	 
	    }
	  }
      }
    /*
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& ibFaces = mesh.getIBFaces();
	  const StorageSite& fIbFaces = fMesh.getIBFaces();
        
	  GeomFields::SSPair key1(&fIbFaces,&fCells);
	  const IMatrix& mIC =
	    dynamic_cast<const IMatrix&>
	    (*_finestGeomFields._interpolationMatrices[key1]);
	      
	  IMatrixV3 mICV3(mIC);  

	  GeomFields::SSPair key2(&fIbFaces,&solidFaces);
	  const IMatrix& mIP =
	    dynamic_cast<const IMatrix&>
	    (*_finestGeomFields._interpolationMatrices[key2]);

	  IMatrixV3 mIPV3(mIP);    

	  shared_ptr<VectorT3Array> ibVvel(new VectorT3Array(ibFaces.getCount()));
	  

	  const VectorT3Array& cVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	  const VectorT3Array& sVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	  
	  ibVvel->zero(); 

	  //velocity interpolation (cells+solidfaces) 
	  mICV3.multiplyAndAdd(*ibVvel,cVel);
	  mIPV3.multiplyAndAdd(*ibVvel,sVel);
	  _macroFields.velocity.addArray(ibFaces,ibVvel);


	}
      }
    */
    }
    
    if (method==2){
      const int numMeshes = _meshes.size();
	const int nSolidFaces = solidFaces.getCount();
 
	shared_ptr<TArray> muSolid(new TArray(nSolidFaces));
	*muSolid =0;
	_macroFields.viscosity.addArray(solidFaces,muSolid);

	shared_ptr<TArray> nueSolid(new TArray(nSolidFaces));
	*nueSolid =0;
	_macroFields.collisionFrequency.addArray(solidFaces,nueSolid);

	const T rho_init=_options["rho_init"]; 
	const T T_init= _options["T_init"]; 
	const T mu_w= _options["mu_w"];
	const T Tmuref= _options["Tmuref"];
	const T muref= _options["muref"];
	const T R=8314.0/_options["molecularWeight"];
	const T nondim_length=_options["nonDimLt"];

	const T mu0=rho_init*R* T_init*nondim_length/pow(2*R* T_init,0.5);  
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[solidFaces]);
	TArray& viscosity = dynamic_cast<TArray&>(_macroFields.viscosity[solidFaces]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[solidFaces]);
	TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[solidFaces]);

	for(int c=0; c<nSolidFaces;c++)
	  {
	    viscosity[c]= muref*pow(temperature[c]*T_init/ Tmuref,mu_w); // viscosity power law
	    collisionFrequency[c]=density[c]*temperature[c]/viscosity[c]*mu0;
	  }
	
	if(_options.fgamma==2){
	  for(int c=0; c<nSolidFaces;c++)
	    collisionFrequency[c]=_options.Prandtl*collisionFrequency[c];
	}
	
 


   for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const Mesh& fMesh = *_finestMeshes[n];
	if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& fCells = fMesh.getCells();
	  const StorageSite& ibFaces = mesh.getIBFaces();
	  const StorageSite& fIbFaces = fMesh.getIBFaces();        

	  GeomFields::SSPair key1(&fIbFaces,&fCells);
	  const IMatrix& mIC =
	    dynamic_cast<const IMatrix&>
	    (*_finestGeomFields._interpolationMatrices[key1]);
	      
	  IMatrix mICV(mIC);
	  IMatrixV3 mICV3(mIC);  

	  GeomFields::SSPair key2(&fIbFaces,&solidFaces);
	  const IMatrix& mIP =
	    dynamic_cast<const IMatrix&>
	    (*_finestGeomFields._interpolationMatrices[key2]);

	  IMatrix mIPV(mIP);
	  IMatrixV3 mIPV3(mIP);    

	  shared_ptr<TArray> ibVtemp(new TArray(ibFaces.getCount()));
	  shared_ptr<TArray> ibVnue(new TArray(ibFaces.getCount()));
	  shared_ptr<TArray> ibVdensity(new TArray(ibFaces.getCount()));
	  shared_ptr<VectorT3Array> ibVvel(new VectorT3Array(ibFaces.getCount()));
	  
	  const TArray& cTemp  = 
	    dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	  const VectorT3Array& cVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	  const TArray& cDensity  = 
	    dynamic_cast<TArray&>(_macroFields.density[cells]);
	  const TArray& sDensity  = 
	    dynamic_cast<TArray&>(_macroFields.density[solidFaces]);
	  const TArray& cNue  = 
	    dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
	  const TArray& sNue  = 
	    dynamic_cast<TArray&>(_macroFields.collisionFrequency[solidFaces]);
	  const TArray& sTemp  = 
	    dynamic_cast<TArray&>(_macroFields.temperature[solidFaces]);
	  const VectorT3Array& sVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	  
	  ibVnue->zero(); 
	  ibVtemp->zero();
	  ibVvel->zero(); 
	  ibVdensity->zero();


	  Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
	  const Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<const Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);

	  TArray cFTemp(fCells.getCount());
	  VectorT3Array cFVel(fCells.getCount());
	  TArray cFDensity(fCells.getCount());
	  TArray cFNue(fCells.getCount());

           for(int c=0;c<fCells.getCount();c++)
	   {
	       cFTemp[c]=cTemp[FinestToCoarse[c][_level]];
	       cFVel[c]=cVel[FinestToCoarse[c][_level]];
	       cFDensity[c]=cDensity[FinestToCoarse[c][_level]];
	       cFNue[c]=cNue[FinestToCoarse[c][_level]];
	   }

	   //nue interpolation (cells)
	  mICV.multiplyAndAdd(*ibVnue,cFNue);
	  mIPV.multiplyAndAdd(*ibVnue,sNue);
	  _macroFields.collisionFrequency.addArray(ibFaces,ibVnue);
	  //temperature interpolation (cells+solidfaces)         
	  mICV.multiplyAndAdd(*ibVtemp,cFTemp);
	  mIPV.multiplyAndAdd(*ibVtemp,sTemp);
	  _macroFields.temperature.addArray(ibFaces,ibVtemp);
	  //density interpolation (cells+solidfaces)         
	  mICV.multiplyAndAdd(*ibVdensity,cFDensity);
	  mIPV.multiplyAndAdd(*ibVdensity,sDensity);
	  _macroFields.density.addArray(ibFaces,ibVdensity);
	  //velocity interpolation (cells+solidfaces) 
	  mICV3.multiplyAndAdd(*ibVvel,cFVel);
	  mIPV3.multiplyAndAdd(*ibVvel,sVel);
	  _macroFields.velocity.addArray(ibFaces,ibVvel);


	}
      }

   const int f_out = 3;
   if (f_out ==1){
     //Step 2 Find fgamma using macroparameters
     const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	    const int numDirections = _quadrature.getDirCount();
	    const StorageSite& ibFaces = mesh.getIBFaces();
	    const int nibFaces=ibFaces.getCount();
	    const double pi=_options.pi;
	    const TArray& ibTemp  =
	      dynamic_cast<TArray&>(_macroFields.temperature[ibFaces]);
	    const VectorT3Array& ibVel =
	      dynamic_cast<VectorT3Array&>(_macroFields.velocity[ibFaces]);
	    const TArray& ibDensity  =
	      dynamic_cast<TArray&>(_macroFields.density[ibFaces]);

	    for (int j=0; j<numDirections; j++)
	      {
		shared_ptr<TArray> ibFndPtrEqES(new TArray(nibFaces));
		TArray&  ibFndEqES= *ibFndPtrEqES;
		ibFndPtrEqES->zero();	
		Field& fndEqES = *_dsfEqPtrES.dsf[j];

		for (int i=0; i<nibFaces; i++)
		  {
		    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
		    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
		    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
		    const T ibu = ibVel[i][0];
		    const T ibv = ibVel[i][1];
		    const T ibw = ibVel[i][2];
		    ibFndEqES[i]=ibDensity[i]/pow(pi*ibTemp[i],1.5)*exp(-(pow(cx[j]-ibu,2.0)+pow(cy[j]-ibv,2.0)+pow(cz[j]-ibw,2.0))/ibTemp[i]);
		  }
		fndEqES.addArray(ibFaces,ibFndPtrEqES);
	      }
	  }
	}
    }
    else if(f_out==2)
      {
	//Step 2 Find fgamma using interpolation (only ESBGK for now)
	for (int n=0; n<numMeshes; n++)
	  {
	    const int numFields= _quadrature.getDirCount();
	    for (int direction = 0; direction < numFields; direction++)
	  {
	    Field& fndEqES = *_dsfEqPtrES.dsf[direction];
	    const Mesh& mesh = *_meshes[n];
	    if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	      const StorageSite& cells = mesh.getCells();
	      const StorageSite& ibFaces = mesh.getIBFaces();
        
	      GeomFields::SSPair key1(&ibFaces,&cells);
	      const IMatrix& mIC =
		dynamic_cast<const IMatrix&>
		(*_geomFields._interpolationMatrices[key1]);
	      
	      IMatrix mICV(mIC);

	      GeomFields::SSPair key2(&ibFaces,&solidFaces);
	      const IMatrix& mIP =
		dynamic_cast<const IMatrix&>
		(*_geomFields._interpolationMatrices[key2]);
	      
	      IMatrix mIPV(mIP);

	      shared_ptr<TArray> ibVf(new TArray(ibFaces.getCount()));
	  
	      const TArray& cf =
		dynamic_cast<const TArray&>(fndEqES[cells]);
	  
	      ibVf->zero();
	      
	      //distribution function interpolation (cells)
	      mICV.multiplyAndAdd(*ibVf,cf);
	      fndEqES.addArray(ibFaces,ibVf);

	    }
	  }
	  }
      }
       //Step3: Relax Distribution function from ibfaces to solid face
    for (int n=0; n<numMeshes; n++)
      {
	const int numDirections = _quadrature.getDirCount();
	for (int j=0; j<numDirections; j++)
	    {
	      Field& fnd = *_dsfPtr.dsf[j];
	      TArray& dsf = dynamic_cast< TArray&>(fnd[solidFaces]);

#ifdef FVM_PARALLEL
	      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE,dsf.getData(),solidFaces.getCount() , MPI::DOUBLE, MPI::SUM);
#endif
	      const Mesh& mesh = *_meshes[n];
	      const Mesh& fMesh = *_finestMeshes[n];
	      if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
		const StorageSite& ibFaces = mesh.getIBFaces();
		const StorageSite& fIbFaces = fMesh.getIBFaces();
		TArray& dsfIB = dynamic_cast< TArray&>(fnd[ibFaces]);
		Field& fndEqES = *_dsfEqPtrES.dsf[j];
		TArray& dsfEqES = dynamic_cast< TArray&>(fndEqES[ibFaces]);
		const StorageSite& faces = fMesh.getFaces();
		const StorageSite& cells = mesh.getCells();
		const CRConnectivity& faceCells = mesh.getAllFaceCells();
		const CRConnectivity& ibFacesTosolidFaces
		  = fMesh.getConnectivity(fIbFaces,solidFaces);
		const IntArray& ibFaceIndices = fMesh.getIBFaceList();
		const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
		const IntArray& sFCRow = ibFacesTosolidFaces.getRow();
		const IntArray& sFCCol = ibFacesTosolidFaces.getCol();
		const int nibFaces = ibFaces.getCount();
		const int nFaces = faces.getCount();
		const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
		const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
		const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    		VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);

		for(int f=0; f<nibFaces; f++)
		  {
		    dsfIB[f]=0.0;				  
		    double distIBSolidInvSum(0.0);
		    for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		      {
			const int c = sFCCol[nc];
			const int faceIB= ibFaceIndices[f];
			const VectorT3Array& solidFaceCentroid =
			  dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[solidFaces]);
			const VectorT3Array& faceCentroid =
			  dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[faces]);
			      
			double distIBSolid (0.0);
			// based on distance - will be thought
			distIBSolid = sqrt(pow((faceCentroid[faceIB][0]-solidFaceCentroid[c][0]),2)+
					   pow((faceCentroid[faceIB][1]-solidFaceCentroid[c][1]),2)+
					   pow((faceCentroid[faceIB][2]-solidFaceCentroid[c][2]),2));
			distIBSolidInvSum += 1/pow(distIBSolid,RelaxDistribution);
		      }
		    for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		      {
			const int c = sFCCol[nc];
			const VectorT3Array& solidFaceCentroid =
			  dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[solidFaces]);
			const VectorT3Array& faceCentroid =
			  dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[faces]);
			const TArray& nue  =
			  dynamic_cast<TArray&>(_macroFields.collisionFrequency[ibFaces]);
			const int faceIB= ibFaceIndices[f];
			// const T coeff = iCoeffs[nc];
			double time_to_wall (0.0);
			double distIBSolid (0.0);
			const T uwall = v[c][0];
			const T vwall = v[c][1];
			const T wwall = v[c][2];
			// based on distance - will be thought
			distIBSolid = sqrt(pow((faceCentroid[faceIB][0]-solidFaceCentroid[c][0]),2)+
					   pow((faceCentroid[faceIB][1]-solidFaceCentroid[c][1]),2)+
					   pow((faceCentroid[faceIB][2]-solidFaceCentroid[c][2]),2));
			time_to_wall = -1*(pow(distIBSolid,2)/((cx[j]-uwall)*(faceCentroid[faceIB][0]-solidFaceCentroid[c][0])+(cy[j]-vwall)*(faceCentroid[faceIB][1]-solidFaceCentroid[c][1])+(cz[j]-wwall)*(faceCentroid[faceIB][2]-solidFaceCentroid[c][2])));
			if(time_to_wall<0)
			  time_to_wall = 0;
	      
			dsfIB[f] += (dsfEqES[f]-(dsfEqES[f]-dsf[c])*exp(-time_to_wall*nue[f]))/(pow(distIBSolid,RelaxDistribution)*distIBSolidInvSum);

		      }

		  }
	      }
	    }
      }
    }
    if (method==3){
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
	      
	  IMatrixV3 mICV3(mIC);  

	  GeomFields::SSPair key2(&ibFaces,&solidFaces);
	  const IMatrix& mIP =
	    dynamic_cast<const IMatrix&>
	    (*_geomFields._interpolationMatrices[key2]);

	  IMatrixV3 mIPV3(mIP);    

	  shared_ptr<VectorT3Array> ibVvel(new VectorT3Array(ibFaces.getCount()));
	  

	  const VectorT3Array& cVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	  const VectorT3Array& sVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	  
	  ibVvel->zero(); 

	  //velocity interpolation (cells+solidfaces) 
	  mICV3.multiplyAndAdd(*ibVvel,cVel);
	  mIPV3.multiplyAndAdd(*ibVvel,sVel);
	  _macroFields.velocity.addArray(ibFaces,ibVvel);


	}
      }
      const int nSolidFaces = solidFaces.getCount();
 
	shared_ptr<TArray> muSolid(new TArray(nSolidFaces));
	*muSolid =0;
	_macroFields.viscosity.addArray(solidFaces,muSolid);

	shared_ptr<TArray> nueSolid(new TArray(nSolidFaces));
	*nueSolid =0;
	_macroFields.collisionFrequency.addArray(solidFaces,nueSolid);

	const T rho_init=_options["rho_init"]; 
	const T T_init= _options["T_init"]; 
	const T mu_w= _options["mu_w"];
	const T Tmuref= _options["Tmuref"];
	const T muref= _options["muref"];
	const T R=8314.0/_options["molecularWeight"];
	const T nondim_length=_options["nonDimLt"];

	const T mu0=rho_init*R* T_init*nondim_length/pow(2*R* T_init,0.5);  
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[solidFaces]);
	TArray& viscosity = dynamic_cast<TArray&>(_macroFields.viscosity[solidFaces]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[solidFaces]);
	TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[solidFaces]);

	for(int c=0; c<nSolidFaces;c++)
	  {
	    viscosity[c]= muref*pow(temperature[c]*T_init/ Tmuref,mu_w); // viscosity power law
	    collisionFrequency[c]=density[c]*temperature[c]/viscosity[c]*mu0;
	  }
	
	if(_options.fgamma==2){
	  for(int c=0; c<nSolidFaces;c++)
	    collisionFrequency[c]=_options.Prandtl*collisionFrequency[c];
	}
	
	//Step 2 Find fgamma using interpolation (only ESBGK for now)
	  const int numFields= _quadrature.getDirCount();
	  for (int direction = 0; direction < numFields; direction++)
	    {
	      shared_ptr<TArray> ibVf(new TArray(solidFaces.getCount()));
	      Field& fndEqES = *_dsfEqPtrES.dsf[direction];
	      ibVf->zero();
	      for (int n=0; n<numMeshes; n++)
		{
		  const Mesh& mesh = *_meshes[n];
		  if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

		    const StorageSite& cells = mesh.getCells();
		    const StorageSite& ibFaces = mesh.getIBFaces();
        
		    GeomFields::SSPair key1(&solidFaces,&cells);
		    const IMatrix& mIC =
		      dynamic_cast<const IMatrix&>
		      (*_geomFields._interpolationMatrices[key1]);
	      
		    IMatrix mICV(mIC);
 
		    const TArray& cf =
		      dynamic_cast<const TArray&>(fndEqES[cells]);
	  
		    ibVf->zero();
	      
		    //distribution function interpolation (cells)
		    mICV.multiplyAndAdd(*ibVf,cf);      
		  }
		}
	      fndEqES.addArray(solidFaces,ibVf);
	    }
    for (int n=0; n<numMeshes; n++)
      {
	const int numDirections = _quadrature.getDirCount();
	for (int j=0; j<numDirections; j++)
	    {
	      Field& fnd = *_dsfPtr.dsf[j];
	      TArray& dsf = dynamic_cast< TArray&>(fnd[solidFaces]); 
	      Field& fndEqES = *_dsfEqPtrES.dsf[j];
	      TArray& dsfEqES = dynamic_cast< TArray&>(fndEqES[solidFaces]);
#ifdef FVM_PARALLEL
	      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE,dsf.getData(),solidFaces.getCount() , MPI::DOUBLE, MPI::SUM);
	      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE,dsfEqES.getData(),solidFaces.getCount() , MPI::DOUBLE, MPI::SUM);
#endif

	      const Mesh& mesh = *_meshes[n];
	      if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
		const StorageSite& ibFaces = mesh.getIBFaces();
		const StorageSite& faces = mesh.getFaces();
		const StorageSite& cells = mesh.getCells();
		const CRConnectivity& faceCells = mesh.getAllFaceCells();
		const CRConnectivity& ibFacesTosolidFaces
		  = mesh.getConnectivity(ibFaces,solidFaces);
		const IntArray& ibFaceIndices = mesh.getIBFaceList();
		const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
		const IntArray& sFCRow = ibFacesTosolidFaces.getRow();
		const IntArray& sFCCol = ibFacesTosolidFaces.getCol();
		const int nibFaces = ibFaces.getCount();
		const int nFaces = faces.getCount();
		const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
		const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
		const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
		VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
		shared_ptr<TArray> ibVf(new TArray(ibFaces.getCount()));
		ibVf->zero();
		TArray&  ibVfA= *ibVf;
		for(int f=0; f<nibFaces; f++)
		  {		  
		    double distIBSolidInvSum(0.0);
		    for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		      {
			const int c = sFCCol[nc];
			const int faceIB= ibFaceIndices[f];
			const VectorT3Array& solidFaceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[solidFaces]);
			const VectorT3Array& faceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[faces]);
			      
			double distIBSolid (0.0);
			// based on distance - will be thought
			distIBSolid = sqrt(pow((faceCentroid[faceIB][0]-solidFaceCentroid[c][0]),2)+
					   pow((faceCentroid[faceIB][1]-solidFaceCentroid[c][1]),2)+
					   pow((faceCentroid[faceIB][2]-solidFaceCentroid[c][2]),2));
			distIBSolidInvSum += 1/pow(distIBSolid,RelaxDistribution);
		      }
		    for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		      {
			const int c = sFCCol[nc];
			const VectorT3Array& solidFaceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[solidFaces]);
			const VectorT3Array& faceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[faces]);
			const TArray& nue  =
			  dynamic_cast<TArray&>(_macroFields.collisionFrequency[solidFaces]);
			const TArray& nueC  =
			  dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
			const int faceIB= ibFaceIndices[f];
			const T uwall = v[c][0];
			const T vwall = v[c][1];
			const T wwall = v[c][2];
			// const T coeff = iCoeffs[nc];
			double time_to_wall (0.0);
			double distIBSolid (0.0);
			// based on distance - will be thought
			distIBSolid = sqrt(pow((faceCentroid[faceIB][0]-solidFaceCentroid[c][0]),2)+
					   pow((faceCentroid[faceIB][1]-solidFaceCentroid[c][1]),2)+
					   pow((faceCentroid[faceIB][2]-solidFaceCentroid[c][2]),2));
			time_to_wall = -1*(pow(distIBSolid,2)/((cx[j]-uwall)*(faceCentroid[faceIB][0]-solidFaceCentroid[c][0])+(cy[j]-vwall)*(faceCentroid[faceIB][1]-solidFaceCentroid[c][1])+(cz[j]-wwall)*(faceCentroid[faceIB][2]-solidFaceCentroid[c][2])));
			if(time_to_wall<0)
			  time_to_wall = 0;
			ibVfA[f] += (dsfEqES[c]-(dsfEqES[c]-dsf[c])*exp(-time_to_wall*nue[c]))/(pow(distIBSolid,RelaxDistribution)*distIBSolidInvSum);
		      }

		  }

		fnd.addArray(ibFaces,ibVf);
	      }
	    }
      }
    }
  }
     
  void computeSolidFacePressure(const StorageSite& solidFaces)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	shared_ptr<TArray> ibP(new TArray(solidFaces.getCount()));
	ibP->zero(); 
	const Mesh& mesh = *_meshes[n];
	if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	  const StorageSite& cells = mesh.getCells();
  
	  GeomFields::SSPair key1(&solidFaces,&cells);
	  const IMatrix& mIC =
	    dynamic_cast<const IMatrix&>
	    (*_geomFields._interpolationMatrices[key1]);
	  
	  IMatrix mICV(mIC);
    
	  const TArray& cP  = 
	    dynamic_cast<TArray&>(_macroFields.pressure[cells]);

	  ibP->zero(); 


	   //nue interpolation (cells)
           mICV.multiplyAndAdd(*ibP,cP);
	}
	_macroFields.pressure.addArray(solidFaces,ibP);
      }
#ifdef FVM_PARALLEL
    TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[solidFaces]);
    MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE,pressure.getData(),solidFaces.getCount() , MPI::DOUBLE, MPI::SUM); 
#endif 
  }

  void computeSolidFaceDsf(const StorageSite& solidFaces,const int method,const int RelaxDistribution=0)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    typedef CRMatrixTranspose<T,VectorT3,VectorT3> IMatrixV3;
    const int numFields= _quadrature.getDirCount();
    if (method==1){
      const int numMeshes = _meshes.size();
      for (int direction = 0; direction < numFields; direction++) {
	Field& fnd = *_dsfPtr.dsf[direction]; 
	shared_ptr<TArray> ibV(new TArray(solidFaces.getCount()));
	ibV->zero();  
	for (int n=0; n<numMeshes; n++)
	  {
	    const Mesh& mesh = *_meshes[n];
	    const Mesh& fMesh = *_finestMeshes[n];
	    if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	    
	      const StorageSite& cells = mesh.getCells();
	      const StorageSite& fCells = fMesh.getCells();
	    
	      GeomFields::SSPair key1(&solidFaces,&fCells);
	      const IMatrix& mIC =
		dynamic_cast<const IMatrix&>
		(*_finestGeomFields._interpolationMatrices[key1]);
	      
	      IMatrix mICV(mIC);

              Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
              const Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<const Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);

	      const TArray& cV =
		dynamic_cast<const TArray&>(fnd[cells]);
              TArray cFV(fCells.getCountLevel1());
              for(int c=0;c<fCells.getCountLevel1();c++)
                cFV[c]=cV[FinestToCoarse[c][_level]];

	      ibV->zero();        
	     
	      mICV.multiplyAndAdd(*ibV,cFV);
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

	    }

	  }
	fnd.addArray(solidFaces,ibV);
      }
    }
    if (method==2){

   // Step0: Compute Interpolation Matrices from (only) Cells to IBFaces
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const Mesh& fMesh = *_finestMeshes[n];       
   
	const int numFields= _quadrature.getDirCount();
	for (int direction = 0; direction < numFields; direction++)
	  {
	    Field& fnd = *_dsfPtr.dsf[direction];
	    Field& fndEqES = *_dsfEqPtrES.dsf[direction];
	    const Mesh& mesh = *_meshes[n];
	    if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	
	      const StorageSite& cells = mesh.getCells();
	      const StorageSite& fCells = fMesh.getCells();
	      const StorageSite& faces = mesh.getFaces();
	      const StorageSite& fFaces = fMesh.getFaces();
	      const StorageSite& ibFaces = mesh.getIBFaces();
	      const StorageSite& fIbFaces = fMesh.getIBFaces();	

	      //GeomFields::SSPair key1(&faces,&cells);
	      GeomFields::SSPair key1(&fFaces,&fCells);
	      //GeomFields::SSPair key1(&fIbFaces,&fCells);
	      const IMatrix& mIC =
		dynamic_cast<const IMatrix&>
		(*_finestGeomFields._interpolationMatrices[key1]);
	      
	      IMatrix mICV(mIC);
	  
	      const TArray& cf =
		dynamic_cast<const TArray&>(fnd[cells]);
	      const TArray& cfEq =
		dynamic_cast<const TArray&>(fndEqES[cells]);

	      Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
	      const Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<const Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);

	      TArray cFf(fCells.getCount());
	      TArray cFfEq(fCells.getCount());

	      for(int c=0;c<fCells.getCount();c++)
	      {
		  cFf[c]=cf[FinestToCoarse[c][_level]];
		  cFfEq[c]=cfEq[FinestToCoarse[c][_level]];
	      }

	      shared_ptr<TArray> ibVf(new TArray(ibFaces.getCount()));
	      ibVf->zero();
	      if (_options.fgamma==2){
		shared_ptr<TArray> ibVfEq(new TArray(ibFaces.getCount()));
		ibVfEq->zero();
		mICV.multiplyAndAdd(*ibVfEq,cFfEq);
		fndEqES.addArray(ibFaces,ibVfEq);
	      }
	      mICV.multiplyAndAdd(*ibVf,cFf);
	      fnd.addArray(ibFaces,ibVf);
           
      }
    }
      }
      const int nSolidFaces = solidFaces.getCount();
 
      shared_ptr<TArray> muSolid(new TArray(nSolidFaces));
      *muSolid =0;
      _macroFields.viscosity.addArray(solidFaces,muSolid);
      
      shared_ptr<TArray> nueSolid(new TArray(nSolidFaces));
      *nueSolid =0;
      _macroFields.collisionFrequency.addArray(solidFaces,nueSolid);

      const T rho_init=_options["rho_init"]; 
      const T T_init= _options["T_init"]; 
      const T mu_w= _options["mu_w"];
      const T Tmuref= _options["Tmuref"];
      const T muref= _options["muref"];
      const T R=8314.0/_options["molecularWeight"];
      const T nondim_length=_options["nonDimLt"];
      
      const T mu0=rho_init*R* T_init*nondim_length/pow(2*R* T_init,0.5);  
	
      TArray& density = dynamic_cast<TArray&>(_macroFields.density[solidFaces]);
      TArray& viscosity = dynamic_cast<TArray&>(_macroFields.viscosity[solidFaces]);
      TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[solidFaces]);
      TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[solidFaces]);

      for(int c=0; c<nSolidFaces;c++)
	{
	  viscosity[c]= muref*pow(temperature[c]*T_init/ Tmuref,mu_w); // viscosity power law
	  collisionFrequency[c]=density[c]*temperature[c]/viscosity[c]*mu0;
	}
	
      if(_options.fgamma==2){
	for(int c=0; c<nSolidFaces;c++)
	  collisionFrequency[c]=_options.Prandtl*collisionFrequency[c];
      }

     //Step 1 Interpolate Macroparameters and f to IBface
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  const Mesh& fMesh = *_finestMeshes[n];
	  if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	    const StorageSite& cells = mesh.getCells();
	    const StorageSite& fCells = fMesh.getCells();
	    const StorageSite& ibFaces = mesh.getIBFaces();
	    const StorageSite& fIbFaces = fMesh.getIBFaces();        

	  GeomFields::SSPair key1(&fIbFaces,&fCells);
	  const IMatrix& mIC =
	    dynamic_cast<const IMatrix&>
	    (*_finestGeomFields._interpolationMatrices[key1]);
	      
	  IMatrix mICV(mIC);
	  IMatrixV3 mICV3(mIC);  

	  GeomFields::SSPair key2(&fIbFaces,&solidFaces);
	  const IMatrix& mIP =
	    dynamic_cast<const IMatrix&>
	    (*_finestGeomFields._interpolationMatrices[key2]);

	  IMatrix mIPV(mIP);
	  IMatrixV3 mIPV3(mIP);    

	  shared_ptr<TArray> ibVtemp(new TArray(ibFaces.getCount()));
	  shared_ptr<TArray> ibVnue(new TArray(ibFaces.getCount()));
	  shared_ptr<TArray> ibVdensity(new TArray(ibFaces.getCount()));
	  shared_ptr<VectorT3Array> ibVvel(new VectorT3Array(ibFaces.getCount()));
	  
	  const TArray& cTemp  = 
	    dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	  const VectorT3Array& cVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	  const TArray& cDensity  = 
	    dynamic_cast<TArray&>(_macroFields.density[cells]);
	  const TArray& sDensity  = 
	    dynamic_cast<TArray&>(_macroFields.density[solidFaces]);
	  const TArray& cNue  = 
	    dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
	  const TArray& sNue  = 
	    dynamic_cast<TArray&>(_macroFields.collisionFrequency[solidFaces]);
	  const TArray& sTemp  = 
	    dynamic_cast<TArray&>(_macroFields.temperature[solidFaces]);
	  const VectorT3Array& sVel = 
	    dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	  
	  ibVnue->zero(); 
	  ibVtemp->zero();
	  ibVvel->zero(); 
	  ibVdensity->zero(); 


          Field& FinestToCoarseField=_finestGeomFields.finestToCoarse;
          const Array<Vector<int,25> >& FinestToCoarse=dynamic_cast<const Array<Vector<int,25> >&>(FinestToCoarseField[fCells]);

          TArray cFTemp(fCells.getCount());
          VectorT3Array cFVel(fCells.getCount());
          TArray cFDensity(fCells.getCount());
          TArray cFNue(fCells.getCount());

	  for(int c=0;c<fCells.getCount();c++)
	  {
	      cFTemp[c]=cTemp[FinestToCoarse[c][_level]];
	      cFVel[c]=cVel[FinestToCoarse[c][_level]];
	      cFDensity[c]=cDensity[FinestToCoarse[c][_level]];
	      cFNue[c]=cNue[FinestToCoarse[c][_level]];
	  }	  

	   //nue interpolation (cells)
           mICV.multiplyAndAdd(*ibVnue,cFNue);
	   mIPV.multiplyAndAdd(*ibVnue,sNue);
           _macroFields.collisionFrequency.addArray(ibFaces,ibVnue);
	   //temperature interpolation (cells+solidfaces)         
	   mICV.multiplyAndAdd(*ibVtemp,cFTemp);
	   mIPV.multiplyAndAdd(*ibVtemp,sTemp);
           _macroFields.temperature.addArray(ibFaces,ibVtemp);
	   //density interpolation (cells+solidfaces)         
	   mICV.multiplyAndAdd(*ibVdensity,cFDensity);
	   mIPV.multiplyAndAdd(*ibVdensity,sDensity);
           _macroFields.density.addArray(ibFaces,ibVdensity);
	   //velocity interpolation (cells+solidfaces) 
	   mICV3.multiplyAndAdd(*ibVvel,cFVel);
	   mIPV3.multiplyAndAdd(*ibVvel,sVel);
           _macroFields.velocity.addArray(ibFaces,ibVvel);


	  }
	}
      if (_options.fgamma==1){
      //Step 2 Find fgamma using macroparameters
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	    const int numDirections = _quadrature.getDirCount();
	    const StorageSite& ibFaces = mesh.getIBFaces();
	    const int nibFaces=ibFaces.getCount();
	    const double pi=_options.pi;
	    const TArray& ibTemp  =
	      dynamic_cast<TArray&>(_macroFields.temperature[ibFaces]);
	    const VectorT3Array& ibVel =
	      dynamic_cast<VectorT3Array&>(_macroFields.velocity[ibFaces]);
	    const TArray& ibDensity  =
	      dynamic_cast<TArray&>(_macroFields.density[ibFaces]);

	    for (int j=0; j<numDirections; j++)
	      {
		shared_ptr<TArray> ibFndPtrEqES(new TArray(nibFaces));
		TArray&  ibFndEqES= *ibFndPtrEqES;
		
		ibFndPtrEqES->zero();

		Field& fndEqES = *_dsfEqPtrES.dsf[j];

		for (int i=0; i<nibFaces; i++)
		  {
		    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
		    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
		    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
		    const T ibu = ibVel[i][0];
		    const T ibv = ibVel[i][1];
		    const T ibw = ibVel[i][2];
		    ibFndEqES[i]=ibDensity[i]/pow(pi*ibTemp[i],1.5)*exp(-(pow(cx[j]-ibu,2.0)+pow(cy[j]-ibv,2.0)+pow(cz[j]-ibw,2.0))/ibTemp[i]);
		  }
		fndEqES.addArray(ibFaces,ibFndPtrEqES);
	      }
	  }
	}
      }

      //Step3: Relax Distribution function from ibfaces to solid face
	const int numDirections = _quadrature.getDirCount();
	for (int j=0; j<numDirections; j++)
	  {
	    const int nSolidFaces = solidFaces.getCount();
	    shared_ptr<TArray> solidFndPtr(new TArray(nSolidFaces));
	    solidFndPtr->zero(); 
	    TArray&  solidFnd= *solidFndPtr;
	    Field& fnd = *_dsfPtr.dsf[j];
	    Field& fndEqES = *_dsfEqPtrES.dsf[j];
	    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	    for (int n=0; n<numMeshes; n++)
	      {
		const Mesh& mesh = *_meshes[n];
		const Mesh& fMesh = *_finestMeshes[n];
		if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
		  const StorageSite& ibFaces = mesh.getIBFaces();
		  const StorageSite& fIbFaces = fMesh.getIBFaces();
		  const CRConnectivity& solidFacesToibFaces
		    = fMesh.getConnectivity(solidFaces,fIbFaces);
		  const IntArray& ibFaceIndices = fMesh.getIBFaceList();
		  const IntArray& sFCRow = solidFacesToibFaces.getRow();
		  const IntArray& sFCCol = solidFacesToibFaces.getCol();
		  TArray& dsf = dynamic_cast< TArray&>(fnd[ibFaces]);  
		  TArray& dsfEqES = dynamic_cast< TArray&>(fndEqES[ibFaces]);
		  VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	    
		  for(int f=0; f<nSolidFaces; f++)
		    {
		      double distIBSolidInvSum(0.0);
		  for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		    {
		      const StorageSite& faces = fMesh.getFaces();
		      const int c = sFCCol[nc];
		      const int faceIB= ibFaceIndices[c];
		      const VectorT3Array& solidFaceCentroid =
			dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[solidFaces]);
		      const VectorT3Array& faceCentroid =
			dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[faces]);
	      
		      double distIBSolid (0.0);
		      // based on distance - will be thought
		      distIBSolid = sqrt(pow((faceCentroid[faceIB][0]-solidFaceCentroid[f][0]),2)+
					 pow((faceCentroid[faceIB][1]-solidFaceCentroid[f][1]),2)+
					 pow((faceCentroid[faceIB][2]-solidFaceCentroid[f][2]),2));
			distIBSolidInvSum += 1/pow(distIBSolid,RelaxDistribution);
		    }
		  for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		    {
		      const int c = sFCCol[nc];
		      const StorageSite& faces = fMesh.getFaces();
		      const VectorT3Array& solidFaceCentroid =
			dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[solidFaces]);
		      const VectorT3Array& faceCentroid =
			dynamic_cast<const VectorT3Array&> (_finestGeomFields.coordinate[faces]);
		      const int faceIB= ibFaceIndices[c];
		      T time_to_wall (0.0);
		      T distIBSolid (0.0);
		      distIBSolid = sqrt(pow((faceCentroid[faceIB][0]-solidFaceCentroid[f][0]),2)+
					 pow((faceCentroid[faceIB][1]-solidFaceCentroid[f][1]),2)+
					 pow((faceCentroid[faceIB][2]-solidFaceCentroid[f][2]),2));
		      // based on distance - will be thought
		      const T uwall = v[f][0];
		      const T vwall = v[f][1];
		      const T wwall = v[f][2];
		      const TArray& nue  =
			dynamic_cast<TArray&>(_macroFields.collisionFrequency[ibFaces]);
		      time_to_wall = (pow(distIBSolid,2)/((cx[j]-uwall)*(faceCentroid[faceIB][0]-solidFaceCentroid[f][0])+(cy[j]-vwall)*(faceCentroid[faceIB][1]-solidFaceCentroid[f][1])+(cz[j]-wwall)*(faceCentroid[faceIB][2]-solidFaceCentroid[f][2])));
		      if(time_to_wall<0)
			time_to_wall = 0;
			  
		      solidFnd[f] += (dsfEqES[c]-(dsfEqES[c]-dsf[c])*exp(-time_to_wall*nue[c]))/(pow(distIBSolid,RelaxDistribution)*distIBSolidInvSum);
		    }

		    }
		}
	      }
	    fnd.addArray(solidFaces,solidFndPtr);
	  }
    }
    if (method==3){
      const int numDirections = _quadrature.getDirCount();
      for (int j=0; j<numDirections; j++)
	{
	  const int nSolidFaces = solidFaces.getCount();
	  shared_ptr<TArray> solidFndPtr(new TArray(nSolidFaces));
	  solidFndPtr->zero(); 
	  TArray&  solidFnd= *solidFndPtr;
	  Field& fnd = *_dsfPtr.dsf[j];
	  Field& fndEqES = *_dsfEqPtrES.dsf[j];
	  const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	  const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	  const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	  const int numMeshes = _meshes.size();
	  for (int n=0; n<numMeshes; n++)
	    {
	      const Mesh& mesh = *_meshes[n];
	      if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
		const StorageSite& cells = mesh.getCells();
		const StorageSite& ibFaces = mesh.getIBFaces();
		const CRConnectivity& solidFacesToCells
		  = mesh.getConnectivity(solidFaces,cells);
		const IntArray& sFCRow = solidFacesToCells.getRow();
		const IntArray& sFCCol = solidFacesToCells.getCol();
		TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);  
		TArray& dsfEqES = dynamic_cast< TArray&>(fndEqES[cells]);
		VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	    
		for(int f=0; f<nSolidFaces; f++)
		  {
		    double distIBSolidInvSum(0.0);
		    for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		      {
			const StorageSite& faces = mesh.getFaces();
			const int c = sFCCol[nc];
			const VectorT3Array& solidFaceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[solidFaces]);
			const VectorT3Array& faceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[cells]);
	      
			double distIBSolid (0.0);
			// based on distance - will be thought
			distIBSolid = sqrt(pow((faceCentroid[c][0]-solidFaceCentroid[f][0]),2)+
					   pow((faceCentroid[c][1]-solidFaceCentroid[f][1]),2)+
					   pow((faceCentroid[c][2]-solidFaceCentroid[f][2]),2));
			distIBSolidInvSum += 1/pow(distIBSolid,RelaxDistribution);
		      }
		    for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
		      {
			const int c = sFCCol[nc];
			const StorageSite& faces = mesh.getFaces();
			const VectorT3Array& solidFaceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[solidFaces]);
			const VectorT3Array& faceCentroid =
			  dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[cells]);
			T time_to_wall (0.0);
			T distIBSolid (0.0);
			const T uwall = v[f][0];
			const T vwall = v[f][1];
			const T wwall = v[f][2];
			distIBSolid = sqrt(pow((faceCentroid[c][0]-solidFaceCentroid[f][0]),2)+
					   pow((faceCentroid[c][1]-solidFaceCentroid[f][1]),2)+
					   pow((faceCentroid[c][2]-solidFaceCentroid[f][2]),2));
			// based on distance - will be thought
			
			const TArray& nue  =
			  dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
			time_to_wall = (pow(distIBSolid,2)/((cx[j]-uwall)*(faceCentroid[c][0]-solidFaceCentroid[f][0])+(cy[j]-vwall)*(faceCentroid[c][1]-solidFaceCentroid[f][1])+(cz[j]-wwall)*(faceCentroid[c][2]-solidFaceCentroid[f][2])));
			if(time_to_wall<0)
			  time_to_wall = 0;
			
			solidFnd[f] += (dsfEqES[c]-(dsfEqES[c]-dsf[c])*exp(-time_to_wall*nue[c]))/(pow(distIBSolid,RelaxDistribution)*distIBSolidInvSum);
		      }
		  
		  }
	      }
	    }	
	  fnd.addArray(solidFaces,solidFndPtr);
	  
	}
    }
  }



  
 void correctMassDeficit()
 {

 const int numMeshes = _meshes.size();
 T netFlux(0.0);
     
  for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
	const StorageSite& ibFaces = mesh.getIBFaces();
        const StorageSite& cells = mesh.getCells();
	const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
        const StorageSite& faces = mesh.getFaces();
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	const int numDirections = _quadrature.getDirCount();
	const IntArray& ibFaceIndex = dynamic_cast<const IntArray&>(_geomFields.ibFaceIndex[faces]);
	const VectorT3Array& faceArea =
       	      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	const TArray& faceAreaMag =
	  dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	const CRConnectivity& faceCells = mesh.getAllFaceCells();
	const int nibFaces = ibFaces.getCount();

	for(int f=0; f<nibFaces; f++)
	  {
	    const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);
	    if (((ibType[c0] == Mesh::IBTYPE_FLUID) && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
		((ibType[c1] == Mesh::IBTYPE_FLUID) && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
	      {
		const int ibFace = ibFaceIndex[f];
		if (ibFace < 0)
		  throw CException("invalid ib face index");
		if (ibType[c0] == Mesh::IBTYPE_FLUID)
		  {
		    const VectorT3 en = faceArea[f]/faceAreaMag[f]; 
		    for (int j=0; j<numDirections; j++)
		      {
			const T c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
			Field& fnd = *_dsfPtr.dsf[j];
			TArray& dsf = dynamic_cast<TArray&>(fnd[ibFaces]); 
		
			netFlux -= dsf[f]*c_dot_en*wts[j]/abs(c_dot_en);
		      }
	      
		  }
		else
		  {
		    const VectorT3 en = faceArea[f]/faceAreaMag[f]; 
		    for (int j=0; j<numDirections; j++)
		      {
			const T c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
			Field& fnd = *_dsfPtr.dsf[j];
			TArray& dsf = dynamic_cast<TArray&>(fnd[ibFaces]); 
		
			netFlux += dsf[f]*c_dot_en*wts[j]/abs(c_dot_en);
		      }
		  }	    
	      }
	  }
    }  

#ifdef FVM_PARALLEL
  MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &netFlux, 1, MPI::DOUBLE, MPI::SUM);
#endif	

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
      const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
			             
      for(int c=0; c<cells.getSelfCount(); c++)
	{
	  if (ibType[c] == Mesh::IBTYPE_FLUID){
	    const int numDirections = _quadrature.getDirCount();
	    T cellMass(0.0);
	    for (int j=0; j<numDirections; j++)
	      {
		Field& fnd = *_dsfPtr.dsf[j];
		TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
		cellMass += wts[j]*dsf[c];
	      }

	    for (int j=0; j<numDirections; j++)
	      {
		Field& fnd = *_dsfPtr.dsf[j];
		TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
		Field& feqES = *_dsfEqPtrES.dsf[j]; //for fgamma_2
		TArray& fgam = dynamic_cast< TArray&>(feqES[cells]);
		fgam[c] = fgam[c]*(1+netFlux*cellVolume[c]/cellMass);
		dsf[c] = dsf[c]*(1+netFlux*cellVolume[c]/cellMass);
	      }
	  }

	}
    }
 }

 void correctMassDeficit2(double n1,double n2)
 {

 const int numMeshes = _meshes.size();
 T netFlux(0.0);
     
 netFlux=n2-n1;

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
      const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
			             
      for(int c=0; c<cells.getSelfCount(); c++)
	{
	  if (ibType[c] == Mesh::IBTYPE_FLUID){
	    const int numDirections = _quadrature.getDirCount();
	    T cellMass(0.0);
	    for (int j=0; j<numDirections; j++)
	      {
		Field& fnd = *_dsfPtr.dsf[j];
		TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
		cellMass += wts[j]*dsf[c];
	      }

	    for (int j=0; j<numDirections; j++)
	      {
		Field& fnd = *_dsfPtr.dsf[j];
		TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
		Field& feqES = *_dsfEqPtrES.dsf[j]; //for fgamma_2
		TArray& fgam = dynamic_cast< TArray&>(feqES[cells]);
		fgam[c] = fgam[c]*(1+netFlux*cellVolume[c]/cellMass);
		dsf[c] = dsf[c]*(1+netFlux*cellVolume[c]/cellMass);
	      }
	  }

	}
    }
 }
 const double ConservationofMassCheck()
 {
  const int numMeshes = _meshes.size();
  T ndens_tot(0.0) ;
     
  for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
	const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
        const StorageSite& faces = mesh.getFaces();
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	const int numDirections = _quadrature.getDirCount();
	const IntArray& ibFaceIndex = dynamic_cast<const IntArray&>(_geomFields.ibFaceIndex[faces]);
	const VectorT3Array& faceArea =
       	      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
	const TArray& faceAreaMag =
	  dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	const CRConnectivity& faceCells = mesh.getAllFaceCells();
	const int nFaces = faces.getCount();

	for(int c=0; c<cells.getCountLevel1(); c++)
	  {
	  if (ibType[c] == Mesh::IBTYPE_FLUID)
	    {
	      for (int j=0; j<numDirections; j++)
		{
		  Field& fnd = *_dsfPtr.dsf[j];
		  TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
		  ndens_tot += wts[j]*dsf[c];
		}
	    }
	  }
    }
    cout << "Hello, I have" << ndens_tot << "number density";
    return ndens_tot;
 }

 void ConservationofMFSolid(const StorageSite& solidFaces) const
 {
    const double pi=_options.pi;
    const double epsilon=_options.epsilon_ES;
    const int nSolidFaces = solidFaces.getCount();
    for (int i=0; i<nSolidFaces; i++)
      {
	const int numDirections = _quadrature.getDirCount();
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const VectorT3Array& solidFaceCentroid =
	  dynamic_cast<const VectorT3Array&>(_finestGeomFields.coordinate[solidFaces]);
	const VectorT3Array& solidFaceArea =
	  dynamic_cast<const VectorT3Array&>(_finestGeomFields.area[solidFaces]);
	const TArray& solidFaceAreaMag =
	  dynamic_cast<const TArray&>(_finestGeomFields.areaMag[solidFaces]);
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_finestMacroFields.velocity[solidFaces]);
	TArray& density  = dynamic_cast<TArray&>(_finestMacroFields.density[solidFaces]);
	TArray& temperature  = dynamic_cast<TArray&>(_finestMacroFields.temperature[solidFaces]);
	const T uwall = v[i][0];
	const T vwall = v[i][1];
	const T wwall = v[i][2];
	const T Twall = temperature[i];

	T Nmr(0.0) ;
	T Dmr(0.0) ;
	T incomFlux(0.0);
	for (int j=0; j<numDirections; j++)
	  {		
	    Field& fnd = *_dsfPtr.dsf[j];
	    TArray& dsf = dynamic_cast< TArray&>(fnd[solidFaces]);
	    const VectorT3 en = solidFaceArea[i]/solidFaceAreaMag[i];
	    const T c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	    const T wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	    const T fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	    
	    if (c_dot_en-wallV_dot_en > 0) //incoming
	      {
		Dmr = Dmr - fwall*wts[j]*(c_dot_en-wallV_dot_en);
		incomFlux=incomFlux-dsf[i]*wts[j]*(c_dot_en-wallV_dot_en);
	      }
	    else
	      {
		Nmr = Nmr + dsf[i]*wts[j]*(c_dot_en-wallV_dot_en);
	      }
	  }
	const T nwall = Nmr/Dmr; // wall number density for initializing Maxwellian
		    
	density[i]=nwall;
		    
	for (int j=0; j<numDirections; j++)
	  {
	    Field& fnd = *_dsfPtr.dsf[j];
	    TArray& dsf = dynamic_cast< TArray&>(fnd[solidFaces]);
	    const VectorT3 en = solidFaceArea[i]/solidFaceAreaMag[i];
	    const T c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	    const T wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	    if (c_dot_en-wallV_dot_en > 0)
	      {
		dsf[i] = nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	      }	    
	    else
	      dsf[i]=dsf[i];
	  }
      }
 }

  void doSweeps(const int sweeps, const int num)
  {
    for(int sweepNo=0;sweepNo<sweeps;sweepNo++)
      smooth(num);
  }

  void doSweeps(const int sweeps, const int num, const StorageSite& solidFaces)
  {
    for(int sweepNo=0;sweepNo<sweeps;sweepNo++)
      smooth(num,solidFaces);
  }

  void smooth(const int num,const StorageSite& solidFaces)
  {
    const int numDir=_quadrature.getDirCount();
    const int numMeshes=_meshes.size();
    for(int msh=0;msh<numMeshes;msh++)
    {
        const Mesh& mesh=*_meshes[msh];
        const BCcellArray& BCArray=*(_BCells[msh]);
        const BCfaceArray& BCfArray=*(_BFaces[msh]);
        const BCcellArray& ZCArray=*(_ZCells[msh]);
        COMETESBGKDiscretizer<T> CDisc(mesh,_geomFields,solidFaces,_macroFields,_quadrature,
                                       _dsfPtr,_dsfPtr1,_dsfPtr2,_dsfEqPtrES,_dsfPtrRes,_dsfPtrFAS,
                                       _options["timeStep"],_options.timeDiscretizationOrder,
                                       _options.transient,_options.underRelaxation,_options["rho_init"], 
                                       _options["T_init"],_options["molecularWeight"],
                                       _bcMap,_faceReflectionArrayMap,BCArray,BCfArray,ZCArray);

        CDisc.setfgFinder();

	//Field::syncLocalVectorFields( _dsfPtr.dsf );
#if 1
        MakeParallel();
#endif  
        CDisc.COMETSolve(1,_level); //forward

	//Field::syncLocalVectorFields( _dsfPtr.dsf );
#if 1
        MakeParallel();
#endif

        //callCOMETBoundaryConditions();
        ComputeCOMETMacroparameters();
        ComputeCollisionfrequency();
        //update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
        if (_options.fgamma==0){initializeMaxwellianEq();}
        else{ EquilibriumDistributionBGK();}
        if (_options.fgamma==2){EquilibriumDistributionESBGK();} 
        computeSolidFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
        ConservationofMFSolid(solidFaces);
        computeIBFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
	          
        CDisc.COMETSolve(-1,_level); //reverse
        if((num==1)||(num==0&&_level==0))
	  {
            //Field::syncLocalVectorFields( _dsfPtr.dsf );
#if 1
            MakeParallel();
#endif
            //callCOMETBoundaryConditions();
            ComputeCOMETMacroparameters();
            ComputeCollisionfrequency();
            //update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
            if (_options.fgamma==0){initializeMaxwellianEq();}
            else{ EquilibriumDistributionBGK();}
            if (_options.fgamma==2){EquilibriumDistributionESBGK();}
	    computeSolidFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
	    ConservationofMFSolid(solidFaces);
	    computeIBFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
	  }
      }
  }

  void smooth(const int num)
  {
    const int numDir=_quadrature.getDirCount();
    const int numMeshes=_meshes.size();
    for(int msh=0;msh<numMeshes;msh++)
    {
        const Mesh& mesh=*_meshes[msh];
        const BCcellArray& BCArray=*(_BCells[msh]);
        const BCfaceArray& BCfArray=*(_BFaces[msh]);
	const BCcellArray& ZCArray=*(_ZCells[msh]);
	shared_ptr<StorageSite> solidFaces(new StorageSite(-1));
	COMETESBGKDiscretizer<T> CDisc(mesh,_geomFields,*solidFaces,_macroFields,_quadrature,
				       _dsfPtr,_dsfPtr1,_dsfPtr2,_dsfEqPtrES,_dsfPtrRes,_dsfPtrFAS,
				       _options["timeStep"],_options.timeDiscretizationOrder,
				       _options.transient,_options.underRelaxation,_options["rho_init"], 
				       _options["T_init"],_options["molecularWeight"],
				       _bcMap,_faceReflectionArrayMap,BCArray,BCfArray,ZCArray);

        CDisc.setfgFinder();
	
	//Field::syncLocalVectorFields( _dsfPtr.dsf );	
#if 1
	MakeParallel();
#endif  
	CDisc.COMETSolve(1,_level); //forward
	//callCOMETBoundaryConditions();
	ComputeCOMETMacroparameters();
	ComputeCollisionfrequency();
	//update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
	if (_options.fgamma==0){initializeMaxwellianEq();}
	else{ EquilibriumDistributionBGK();}
	if (_options.fgamma==2){EquilibriumDistributionESBGK();} 
      
        //Field::syncLocalVectorFields( _dsfPtr.dsf );
#if 1        
	MakeParallel();
#endif          
	CDisc.COMETSolve(-1,_level); //reverse
	if((num==1)||(num==0&&_level==0))
	{
	    //callCOMETBoundaryConditions();
	    ComputeCOMETMacroparameters();
	    ComputeCollisionfrequency();
	    //update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
	    if (_options.fgamma==0){initializeMaxwellianEq();}
	    else{ EquilibriumDistributionBGK();}
	    if (_options.fgamma==2){EquilibriumDistributionESBGK();}
	}
      }
  }

  T updateResid(const bool addFAS)
  {
    const int numMeshes=_meshes.size();
    T lowResid=-1.;
    T currentResid;
    for(int msh=0;msh<numMeshes;msh++)
    {
        const Mesh& mesh=*_meshes[msh];
        const BCcellArray& BCArray=*(_BCells[msh]);
        const BCfaceArray& BCfArray=*(_BFaces[msh]);
        const BCcellArray& ZCArray=*(_ZCells[msh]);
	shared_ptr<StorageSite> solidFaces(new StorageSite(-1));
        COMETESBGKDiscretizer<T> CDisc(mesh,_geomFields,*solidFaces,_macroFields,_quadrature,
                                       _dsfPtr,_dsfPtr1,_dsfPtr2,_dsfEqPtrES,_dsfPtrRes,_dsfPtrFAS,
                                       _options["timeStep"],_options.timeDiscretizationOrder,
                                       _options.transient,_options.underRelaxation,_options["rho_init"], 
                                       _options["T_init"],_options["molecularWeight"],
                                       _bcMap,_faceReflectionArrayMap,BCArray,BCfArray,ZCArray);

        CDisc.setfgFinder();
        const int numDir=_quadrature.getDirCount();

	//Field::syncLocalVectorFields( _dsfPtr.dsf );
#if 1
        MakeParallel();
#endif
        CDisc.findResid(addFAS);
        currentResid=CDisc.getAveResid();

        if(lowResid<0)
          lowResid=currentResid;
        else
          if(currentResid<lowResid)
            lowResid=currentResid;
    }
    return lowResid;
  }

  T updateResid(const bool addFAS,const StorageSite& solidFaces)
  {
    const int numMeshes=_meshes.size();
    T lowResid=-1.;
    T currentResid;
    for(int msh=0;msh<numMeshes;msh++)
      {
        const Mesh& mesh=*_meshes[msh];
        const BCcellArray& BCArray=*(_BCells[msh]);
        const BCfaceArray& BCfArray=*(_BFaces[msh]);
	const BCcellArray& ZCArray=*(_ZCells[msh]);
	COMETESBGKDiscretizer<T> CDisc(mesh,_geomFields,solidFaces,_macroFields,_quadrature,
				       _dsfPtr,_dsfPtr1,_dsfPtr2,_dsfEqPtrES,_dsfPtrRes,_dsfPtrFAS,
				       _options["timeStep"],_options.timeDiscretizationOrder,
				       _options.transient,_options.underRelaxation,_options["rho_init"], 
				       _options["T_init"],_options["molecularWeight"],
				       _bcMap,_faceReflectionArrayMap,BCArray,BCfArray,ZCArray);

        CDisc.setfgFinder();
	const int numDir=_quadrature.getDirCount();
	
	//Field::syncLocalVectorFields( _dsfPtr.dsf );
#if 1
	MakeParallel();
#endif
        CDisc.findResid(addFAS);
        currentResid=CDisc.getAveResid();

        if(lowResid<0)
          lowResid=currentResid;
        else
          if(currentResid<lowResid)
            lowResid=currentResid;
      }
    return lowResid;
  }

  void cycle()
  {
    if(_level+1<_options.maxLevels)
      doSweeps(_options.preSweeps,1);
    else
      doSweeps(_options.preCoarsestSweeps,1);
    if(_level+1<_options.maxLevels)
    {
        if(_level==0)
          updateResid(false);
        else
          updateResid(true);

        injectResid();
        _coarserLevel->ComputeCOMETMacroparameters();
        _coarserLevel->ComputeCollisionfrequency();
        if (_options.fgamma==0){_coarserLevel->initializeMaxwellianEq();}
        else{_coarserLevel->EquilibriumDistributionBGK();}
        if (_options.fgamma==2){_coarserLevel->EquilibriumDistributionESBGK();}

        _coarserLevel->makeFAS();
        _coarserLevel->cycle();
        correctSolution();

        ComputeCOMETMacroparameters();
        ComputeCollisionfrequency();
        if (_options.fgamma==0){initializeMaxwellianEq();}
        else{EquilibriumDistributionBGK();}
        if (_options.fgamma==2){EquilibriumDistributionESBGK();}
    }
    if(_level+1<_options.maxLevels)
      doSweeps(_options.postSweeps,0);
    else
    {
	if(_options.postCoarsestSweeps==1)
	  doSweeps(_options.postCoarsestSweeps,0);
        else if(_options.postCoarsestSweeps>1)
          doSweeps(_options.postCoarsestSweeps,1);
    }
  }


  void cycle(const StorageSite& solidFaces)
  {
    if(_level+1<_options.maxLevels)
      doSweeps(_options.preSweeps,1,solidFaces);
    else
      doSweeps(_options.preCoarsestSweeps,1,solidFaces);
    if(_level+1<_options.maxLevels)
      {
        if(_level==0)
          updateResid(false,solidFaces);
        else
          updateResid(true,solidFaces);

        injectResid();
	_coarserLevel->ComputeCOMETMacroparameters();
	_coarserLevel->ComputeCollisionfrequency();
        if (_options.fgamma==0){_coarserLevel->initializeMaxwellianEq();}
	else{_coarserLevel->EquilibriumDistributionBGK();}
	if (_options.fgamma==2){_coarserLevel->EquilibriumDistributionESBGK();}
	_coarserLevel->MakeParallel();     
        _coarserLevel->computeSolidFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
        _coarserLevel->ConservationofMFSolid(solidFaces);
        _coarserLevel->computeIBFaceDsf(solidFaces,_options.method,_options.relaxDistribution);

        _coarserLevel->makeFAS(solidFaces);
        _coarserLevel->cycle(solidFaces);
        correctSolution();
	
        ComputeCOMETMacroparameters();
        ComputeCollisionfrequency();
        if (_options.fgamma==0){initializeMaxwellianEq();}
        else{EquilibriumDistributionBGK();}
        if (_options.fgamma==2){EquilibriumDistributionESBGK();}
	MakeParallel();
        computeSolidFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
        ConservationofMFSolid(solidFaces);
        computeIBFaceDsf(solidFaces,_options.method,_options.relaxDistribution);
      }
    if(_level+1<_options.maxLevels)
      doSweeps(_options.postSweeps,0,solidFaces);
    else
    {
        if(_options.postCoarsestSweeps==1)
          doSweeps(_options.postCoarsestSweeps,0,solidFaces);
        else if(_options.postCoarsestSweeps>1)
          doSweeps(_options.postCoarsestSweeps,1,solidFaces);
    }
  }

  void injectResid()
  {
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& finerMesh=*_meshes[n];
        const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
        MacroFields& coarserMacro=_coarserLevel->getMacro();
	GeomFields& coarserGeomFields=_coarserLevel->getGeomFields();
        DistFunctFields<T>& coarserdsf = _coarserLevel->getdsf();
        DistFunctFields<T>& coarserdsf0 = _coarserLevel->getdsf0();
	DistFunctFields<T>& coarserdsfFAS = _coarserLevel->getdsfFAS();
        const StorageSite& finerCells=finerMesh.getCells();
        const StorageSite& coarserCells=coarserMesh.getCells();
        const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);
        //const TArray& coarserVol=dynamic_cast<TArray&>(_geomFields.volume[coarserCells]);
	const TArray& coarserVol=dynamic_cast<TArray&>(coarserGeomFields.volume[coarserCells]);
        const TArray& finerVol=dynamic_cast<TArray&>(_geomFields.volume[finerCells]);

        const int cellCount=coarserCells.getSelfCount();
        const int numDir=_quadrature.getDirCount();
        for(int dir=0;dir<numDir;dir++)
        {
	    Field& fnd = *_dsfPtr.dsf[dir];
	    Field& fndRes = *_dsfPtrRes.dsf[dir];
	    Field& cfnd = *coarserdsf.dsf[dir];
	    Field& cfndInj = *coarserdsf0.dsf[dir];
	    Field& cfndFAS = *coarserdsfFAS.dsf[dir];
	    TArray& coarserVar=dynamic_cast<TArray&>(cfnd[coarserCells]);
	    TArray& coarserInj=dynamic_cast<TArray&>(cfndInj[coarserCells]);
	    TArray& coarserFAS=dynamic_cast<TArray&>(cfndFAS[coarserCells]);
	    TArray& finerVar=dynamic_cast<TArray&>(fnd[finerCells]);
	    TArray& finerRes=dynamic_cast<TArray&>(fndRes[finerCells]);

 	    for(int c=0;c<cellCount;c++)
	    {
		const int fineCount=CoarserToFiner.getCount(c);
		coarserVar[c]=0.;
		coarserFAS[c]=0.;
                
		for(int fc=0;fc<fineCount;fc++)
		{
		    const int cell=CoarserToFiner(c,fc);
		    coarserVar[c]+=finerVar[cell]*finerVol[cell];
		    coarserFAS[c]+=finerRes[cell];
		}
		coarserVar[c]/=coarserVol[c];
		coarserInj[c]=coarserVar[c];
	    }              
	}

        VectorT3Array& coarserVar=dynamic_cast<VectorT3Array&>(coarserMacro.velocity[coarserCells]);
        VectorT3Array& coarserInj=dynamic_cast<VectorT3Array&>(coarserMacro.velocityInjected[coarserCells]);
        VectorT3Array& coarserFAS=dynamic_cast<VectorT3Array&>(coarserMacro.velocityFASCorrection[coarserCells]);
        VectorT3Array& finerVar=dynamic_cast<VectorT3Array&>(_macroFields.velocity[finerCells]);
        VectorT3Array& finerRes=dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[finerCells]);

        for(int c=0;c<cellCount;c++)
        {
            const int fineCount=CoarserToFiner.getCount(c);
            coarserVar[c]=0.;
            coarserFAS[c]=0.;

            for(int fc=0;fc<fineCount;fc++)
            {
                const int cell=CoarserToFiner(c,fc);
                coarserVar[c]+=finerVar[cell]*finerVol[cell];
                coarserFAS[c]+=finerRes[cell];
	    }
            coarserVar[c]/=coarserVol[c];
            coarserInj[c]=coarserVar[c];
	}
    }
  }

  void makeFAS()
  {
    updateResid(false);

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh=*_meshes[n];
        const StorageSite& cells=mesh.getCells();

        const int numDir = _quadrature.getDirCount();
        for(int dir=0;dir<numDir;dir++)
	{
            Field& fndRes = *_dsfPtrRes.dsf[dir];
            Field& fndFAS = *_dsfPtrFAS.dsf[dir];
            TArray& fRes = dynamic_cast<TArray&>(fndRes[cells]);
            TArray& fFAS = dynamic_cast<TArray&>(fndFAS[cells]);

            fFAS-=fRes;
	}

        VectorT3Array& vR = dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[cells]);
        VectorT3Array& vF = dynamic_cast<VectorT3Array&>(_macroFields.velocityFASCorrection[cells]);

        vF-=vR;
    }
  }

  void makeFAS(const StorageSite& solidFaces)
  {
    updateResid(false,solidFaces);

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh=*_meshes[n];
        const StorageSite& cells=mesh.getCells();

	const int numDir = _quadrature.getDirCount();
        for(int dir=0;dir<numDir;dir++)
        {
	    Field& fndRes = *_dsfPtrRes.dsf[dir];
	    Field& fndFAS = *_dsfPtrFAS.dsf[dir];
	    TArray& fRes = dynamic_cast<TArray&>(fndRes[cells]); 
	    TArray& fFAS = dynamic_cast<TArray&>(fndFAS[cells]);
              
	    fFAS-=fRes;              
	}

        VectorT3Array& vR = dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[cells]);
        VectorT3Array& vF = dynamic_cast<VectorT3Array&>(_macroFields.velocityFASCorrection[cells]);

        vF-=vR;
      }
  }

  void correctSolution()
  {
    const int numMeshes = _meshes.size();
    
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& finerMesh=*_meshes[n];
        const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
        MacroFields& coarserMacro=_coarserLevel->getMacro();
	DistFunctFields<T>& coarserdsf = _coarserLevel->getdsf();
	DistFunctFields<T>& coarserdsf0 = _coarserLevel->getdsf0();
        const StorageSite& finerCells=finerMesh.getCells();
        const StorageSite& coarserCells=coarserMesh.getCells();
        const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);

        const int cellCount=coarserCells.getSelfCount();

        const int numDir=_quadrature.getDirCount();
        for(int dir=0;dir<numDir;dir++)
        {
            Field& fnd = *_dsfPtr.dsf[dir];
	    Field& cfnd = *coarserdsf.dsf[dir];
	    Field& cfndInj = *coarserdsf0.dsf[dir];
            TArray& finerArray = dynamic_cast<TArray&>(fnd[finerCells]);
	    TArray& coarserArray = dynamic_cast<TArray&>(cfnd[coarserCells]);
	    TArray& injArray = dynamic_cast<TArray&>(cfndInj[coarserCells]);
            
	    for(int c=0;c<cellCount;c++)
	    {
		const int fineCount=CoarserToFiner.getCount(c);
		const T correction=coarserArray[c]-injArray[c];
                
		for(int fc=0;fc<fineCount;fc++)
		  finerArray[CoarserToFiner(c,fc)]+=correction;
	    }
	}
	
        VectorT3Array& coarserArray=dynamic_cast<VectorT3Array&>(coarserMacro.velocity[coarserCells]);
        VectorT3Array& injArray=dynamic_cast<VectorT3Array&>(coarserMacro.velocityInjected[coarserCells]);
        VectorT3Array& finerArray=dynamic_cast<VectorT3Array&>(_macroFields.velocity[finerCells]);

        for(int c=0;c<cellCount;c++)
        {
            const int fineCount=CoarserToFiner.getCount(c);
            const VectorT3 correction=coarserArray[c]-injArray[c];
            
            for(int fc=0;fc<fineCount;fc++)
              finerArray[CoarserToFiner(c,fc)]+=correction;
	}
    }
  }

  void advance(const int iters)
  {
    callCOMETBoundaryConditions();
    _residual=updateResid(false);
    _initialResidual=_residual;
    T residualRatio(1.0);
    
#ifdef FVM_PARALLEL    
    if ( MPI::COMM_WORLD.Get_rank() == 0 )
      cout<<"Initial Residual:"<<_initialResidual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif


#ifndef FVM_PARALLEL    
    cout << "Initial Residual:"<<_initialResidual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif
    
    
    
    int niters=0;
    const T absTol=_options.absoluteTolerance;
    const T relTol=_options.relativeTolerance;
    const int show=_options.showResidual;

    while((niters<iters) && ((_residual>absTol)&&(residualRatio>relTol)))
      {
        cycle();
        niters++;
        _residual=updateResid(false);
        if(niters==1)
          _initialResidual=_residual;
        residualRatio=_residual/_initialResidual;
#ifdef FVM_PARALLEL
        if((niters%show==0)&&(MPI::COMM_WORLD.Get_rank()==0))
          cout<<"Iteration:"<<niters<<" Residual:"<<_residual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif

#ifndef FVM_PARALLEL
        if(niters%show==0)
          cout<<"Iteration:"<<niters<<" Residual:"<<_residual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif
      }
    callCOMETBoundaryConditions();
    //cout<<endl<<"Total Iterations:"<<niters<<" Residual:"<<_residual<<endl;
  }

  T advance(const int iters,const StorageSite& solidFaces)
  {
    callCOMETBoundaryConditions();
    _residual=updateResid(false,solidFaces);
    _initialResidual=_residual;
    T residualRatio(1.0);
    
#ifdef FVM_PARALLEL    
    if ( MPI::COMM_WORLD.Get_rank() == 0 )
        cout<<"Initial Residual:"<<_initialResidual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif


#ifndef FVM_PARALLEL    
    cout << "Initial Residual:"<<_initialResidual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif	
    
    
    
    int niters=0;
    const T absTol=_options.absoluteTolerance;
    const T relTol=_options.relativeTolerance;
    const int show=_options.showResidual;

    while((niters<iters) && ((_residual>absTol)&&(residualRatio>relTol)))
      {
        cycle(solidFaces);
        niters++;
        _residual=updateResid(false,solidFaces);
        if(niters==1)
          _initialResidual=_residual;
	residualRatio=_residual/_initialResidual;
#ifdef FVM_PARALLEL
        if((niters%show==0)&&(MPI::COMM_WORLD.Get_rank()==0))
          cout<<"Iteration:"<<niters<<" Residual:"<<_residual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif

#ifndef FVM_PARALLEL
        if(niters%show==0)
          cout<<"Iteration:"<<niters<<" Residual:"<<_residual<<"  ResidualRatio: "<<residualRatio<<endl;
#endif
      }
    callCOMETBoundaryConditions();
    return _residual;
    //cout<<endl<<"Total Iterations:"<<niters<<" Residual:"<<_residual<<endl;
  }
  
  void  computeSurfaceForce(const StorageSite& solidFaces, bool perUnitArea, bool IBM=0)
  {
    typedef CRMatrixTranspose<T,T,T> IMatrix;
    
    const int nSolidFaces = solidFaces.getCount();
    
    boost::shared_ptr<VectorT3Array>
      forcePtr( new VectorT3Array(nSolidFaces));
    VectorT3Array& force = *forcePtr;
    
    force.zero();
    _macroFields.force.addArray(solidFaces,forcePtr);
    
    const VectorT3Array& solidFaceArea =
      dynamic_cast<const VectorT3Array&>(_geomFields.area[solidFaces]);
    
    const TArray& solidFaceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[solidFaces]);
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
        const StorageSite& cells = mesh.getCells();
        
        const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
        const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
        const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
        const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
        const TArray& wts = dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

        const CRConnectivity& solidFacesToCells
          = mesh.getConnectivity(solidFaces,cells);
        const IntArray& sFCRow = solidFacesToCells.getRow();
        const IntArray& sFCCol = solidFacesToCells.getCol(); 

        const T Lx=_options["nonDimLx"];
        const T Ly=_options["nonDimLy"];
        const T Lz=_options["nonDimLz"];
    

        const int N123= _quadrature.getDirCount(); 

        const int selfCount = cells.getSelfCount();
        for(int f=0; f<nSolidFaces; f++){
          
          StressTensor<T> stress = NumTypeTraits<StressTensor<T> >::getZero();
          
          if (IBM){
	    const VectorT3Array& vs = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[solidFaces]);
	    GeomFields::SSPair key1(&solidFaces,&cells);
            const IMatrix& mIC =
              dynamic_cast<const IMatrix&>
              (*_geomFields._interpolationMatrices[key1]);
            const Array<T>& iCoeffs = mIC.getCoeff();
            for(int j=0;j<N123;j++){
              Field& fnd = *_dsfPtr.dsf[j];
              const TArray& f_dsf = dynamic_cast<const TArray&>(fnd[cells]);
              const TArray& fs_dsf = dynamic_cast<const TArray&>(fnd[solidFaces]);
              stress[0] -=pow((cx[j]-vs[f][0]),2.0)*fs_dsf[f]*wts[j];
              stress[1] -=pow((cy[j]-vs[f][1]),2.0)*fs_dsf[f]*wts[j];
              stress[2] -=pow((cz[j]-vs[f][2]),2.0)*fs_dsf[f]*wts[j];
              stress[3] -=(cx[j]-vs[f][0])*(cy[j]-vs[f][1])*fs_dsf[f]*wts[j];
              stress[4] -=(cy[j]-vs[f][1])*(cz[j]-vs[f][2])*fs_dsf[f]*wts[j];
              stress[5] -=(cx[j]-vs[f][0])*(cz[j]-vs[f][2])*fs_dsf[f]*wts[j];
	      /*
              for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
                {
                  
                  const int c = sFCCol[nc];            
                  const T coeff = iCoeffs[nc];
                  stress[0] -=coeff*pow((cx[j]-v[c][0]),2.0)*f_dsf[c]*wts[j];
                  stress[1] -=coeff*pow((cy[j]-v[c][1]),2.0)*f_dsf[c]*wts[j];
                  stress[2] -=coeff*pow((cz[j]-v[c][2]),2.0)*f_dsf[c]*wts[j];
                  stress[3] -=coeff*(cx[j]-v[c][0])*(cy[j]-v[c][1])*f_dsf[c]*wts[j];
                  stress[4] -=coeff*(cy[j]-v[c][1])*(cz[j]-v[c][2])*f_dsf[c]*wts[j];
                  stress[5] -=coeff*(cx[j]-v[c][0])*(cz[j]-v[c][2])*f_dsf[c]*wts[j];
                }
	      */
            }
          }
          else
            {
              for(int j=0;j<N123;j++){
                Field& fnd = *_dsfPtr.dsf[j];
                const TArray& f_dsf = dynamic_cast<const TArray&>(fnd[cells]);
                for(int nc = sFCRow[f]; nc<sFCRow[f+1]; nc++)
                  {
                  
                    const int c = sFCCol[nc];            
                    if ( c < selfCount ){
                      stress[0] +=pow((cx[j]-v[c][0]),2.0)*f_dsf[c]*wts[j];
                      stress[1] +=pow((cy[j]-v[c][1]),2.0)*f_dsf[c]*wts[j];
                      stress[2] +=pow((cz[j]-v[c][2]),2.0)*f_dsf[c]*wts[j];
                      stress[3] +=(cx[j]-v[c][0])*(cy[j]-v[c][1])*f_dsf[c]*wts[j];
                      stress[4] +=(cy[j]-v[c][1])*(cz[j]-v[c][2])*f_dsf[c]*wts[j];
                      stress[5] +=(cx[j]-v[c][0])*(cz[j]-v[c][2])*f_dsf[c]*wts[j];
                    }
                  }
              }

            }

          
          const VectorT3& Af = solidFaceArea[f];
          force[f][0] = Af[0]*Ly*Lz*stress[0] + Af[1]*Lz*Lx*stress[3] + Af[2]*Lx*Ly*stress[5];
          force[f][1] = Af[0]*Ly*Lz*stress[3] + Af[1]*Lz*Lx*stress[1] + Af[2]*Lx*Ly*stress[4];
          force[f][2] = Af[0]*Ly*Lz*stress[5] + Af[1]*Lz*Lx*stress[4] + Af[2]*Ly*Ly*stress[2];
          if (perUnitArea){
            force[f] /= solidFaceAreaMag[f];}
        }
      }
  }
  
  void OutputDsfPOINT() //, const char* filename)
  {
    FILE * pFile;
    pFile = fopen("cxyz0.plt","w");  
    int N1=_quadrature.getNVCount();
    int N2=_quadrature.getNthetaCount();
    int N3=_quadrature.getNphiCount();
    fprintf(pFile,"%s \n", "VARIABLES= cx, cy, cz, f,");
    fprintf(pFile, "%s %i %s %i %s %i \n","ZONE I=", N3,",J=",N2,"K=",N1);
    fprintf(pFile,"%s\n","F=POINT");
    const int numMeshes = _meshes.size();
    const int cellno=_options.printCellNumber;
    for (int n=0; n<numMeshes; n++)
      {
      const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
      const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
      const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
      const int numFields= _quadrature.getDirCount(); 	
      
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      for(int j=0;j< numFields;j++)
	{
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& f = dynamic_cast< TArray&>(fnd[cells]);	
	Field& fEqnd = *_dsfEqPtr.dsf[j];
	TArray& fEq = dynamic_cast< TArray&>(fEqnd[cells]);
	fprintf(pFile,"%E %E %E %E %E\n",cx[j],cy[j],cz[j],f[cellno],fEq[cellno]);
	}
      }
  }
 
  DistFunctFields<T>& getdsf() { return _dsfPtr;} 
  const DistFunctFields<T>& getdsf1() const { return _dsfPtr1;} 
  const DistFunctFields<T>& getdsf2() const { return _dsfPtr2;}
  const DistFunctFields<T>& getdsfEq() const { return _dsfEqPtr;} 
  const DistFunctFields<T>& getdsfEqES() const { return _dsfEqPtrES;}
  DistFunctFields<T>& getdsf0() { return _dsfPtr0;}
  DistFunctFields<T>& getdsfInj() { return _dsfPtrInj;}
  DistFunctFields<T>& getdsfRes() { return _dsfPtrRes;}
  DistFunctFields<T>& getdsfFAS() { return _dsfPtrFAS;}
  void setBCMap(COMETBCMap* bcMap) {_bcMap=*bcMap;}
  void setCoarserLevel(TCOMET* cl) {_coarserLevel=cl;}
  void setFinerLevel(TCOMET* fl) {_finerLevel=fl;}
  int getLevel() {return _level;}
  const MeshList& getMeshList() {return _meshes;}
  GeomFields& getGeomFields() {return _geomFields;}
  TQuad& getQuadrature() {return _quadrature;}
  MacroFields& getMacro() {return _macroFields;}
  T getResidual() {return _residual;}
  
    

 private:
  //shared_ptr<Impl> _impl;
 
  const int _level;
  GeomFields& _geomFields;
  GeomFields _coarseGeomFields;
  GeomFields& _finestGeomFields;
  const MeshList& _finestMeshes;
  Quadrature<T>& _quadrature;
  const int _ibm;
 
  MacroFields& _macroFields;
  MacroFields& _finestMacroFields;
  TCOMET* _finestLevel;
  TCOMET* _coarserLevel;
  TCOMET* _finerLevel;
  DistFunctFields<T> _dsfPtr;  
  DistFunctFields<T> _dsfPtr1;
  DistFunctFields<T> _dsfPtr2;
  DistFunctFields<T> _dsfEqPtr;
  DistFunctFields<T> _dsfEqPtrES;
  DistFunctFields<T> _dsfPtr0;
  DistFunctFields<T> _dsfPtrInj;
  DistFunctFields<T> _dsfPtrRes;
  DistFunctFields<T> _dsfPtrFAS;

  COMETBCMap _bcMap;
  COMETVCMap _vcMap;

  COMETModelOptions<T> _options;
  T _residual;
  T _initialResidual;
  

  MFRPtr _initialKmodelNorm;
  BCcellList _BCells;
  BCfaceList _BFaces;
  BCcellList _ZCells;
  int _niters;
  map<int, vector<int> > _faceReflectionArrayMap;  
  map<string,shared_ptr<ArrayBase> > _persistenceData;

  //MatrixSizeMap _coarseSizes;
  //MatrixSizeMap _coarseGhostSizes;
  SizeMap _coarseSizes;
  SizeMap _coarseGhostSizes;
  SiteMap _siteMap;
  GhostStorageSiteMap _sharedSiteMap;
  MatrixMappersMap _coarseScatterMaps;
  MatrixMappersMap _coarseGatherMaps;
};

#endif
