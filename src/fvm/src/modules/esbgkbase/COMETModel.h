#ifndef _COMETMODEL_H_
#define _COMETMODEL_H_

#include <stdio.h>
#include <map>
#include <cmath>
#include <vector>
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

  /**
   * Calculation of macro-parameters density, temperature, components of velocity, pressure
   * by taking moments of distribution function using quadrature points and weights from quadrature.h
   */
  //MacroFields& macroFields;
  
  COMETModel(const MeshList& meshes, const int level, GeomFields& geomFields, MacroFields& macroFields, Quadrature<T>& quad):
       
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
    _initialKmodelNorm(),
    _niters(0),
    _residual(0.0),
    _initialResidual(0.0)
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

	  const COMETVC<T>& vc = *_vcMap[mesh.getID()];
	  
	  const StorageSite& cells = mesh.getCells();
	  
	  const int nCells = cells.getCount();
	  shared_ptr<VectorT3Array> vCell(new VectorT3Array(nCells));
	  
	  VectorT3 initialVelocity;
	  initialVelocity[0] = _options["initialXVelocity"];
	  initialVelocity[1] = _options["initialYVelocity"];
	  initialVelocity[2] = _options["initialZVelocity"];
	  *vCell = initialVelocity;
	  _macroFields.velocity.addArray(cells,vCell);
	  
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

	  shared_ptr<TArray> tempCell(new TArray(cells.getCount()));
	  *tempCell = _options["operatingTemperature"];
	  _macroFields.temperature.addArray(cells,tempCell);

	  shared_ptr<TArray> collFreqCell(new TArray(cells.getCount()));
	  *collFreqCell = 0.;
	  //*collFreqCell = vc["viscosity"];
	  _macroFields.collisionFrequency.addArray(cells,collFreqCell);
	  
	  //coeffs for perturbed BGK distribution function
	  shared_ptr<VectorT5Array> coeffCell(new VectorT5Array(cells.getCount()));
	  VectorT5 initialCoeff;
	  initialCoeff[0] = 1.0;
	  initialCoeff[1] = 1.0;
	  initialCoeff[2] = 0.0; 
	  initialCoeff[3] = 0.0;
	  initialCoeff[4] = 0.0;
	  *coeffCell = initialCoeff;
	  _macroFields.coeff.addArray(cells,coeffCell);
	  
	  //coeffs for perturbed BGK distribution function
	  shared_ptr<VectorT10Array> coeffgCell(new VectorT10Array(cells.getCount()));
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
	  shared_ptr<TArray> tempxxCell(new TArray(cells.getCount()));
	  *tempxxCell = _options["operatingTemperature"]/3;
	  _macroFields.Txx.addArray(cells,tempxxCell);
	  
	  shared_ptr<TArray> tempyyCell(new TArray(cells.getCount()));
	  *tempyyCell = _options["operatingTemperature"]/3;
	  _macroFields.Tyy.addArray(cells,tempyyCell);
	  
	  shared_ptr<TArray> tempzzCell(new TArray(cells.getCount()));
	  *tempzzCell = _options["operatingTemperature"]/3;
	  _macroFields.Tzz.addArray(cells,tempzzCell);
	  
	  shared_ptr<TArray> tempxyCell(new TArray(cells.getCount()));
	  *tempxyCell = 0.0;
	  _macroFields.Txy.addArray(cells,tempxyCell);
	  
	  shared_ptr<TArray> tempyzCell(new TArray(cells.getCount()));
	  *tempyzCell = 0.0;
	  _macroFields.Tyz.addArray(cells,tempyzCell);
	  
	  shared_ptr<TArray> tempzxCell(new TArray(cells.getCount()));
	  *tempzxCell = 0.0;
	  _macroFields.Tzx.addArray(cells,tempzxCell);
	  
	//Entropy and Entropy Generation Rate for switching
	  shared_ptr<TArray> EntropyCell(new TArray(cells.getCount()));
	  *EntropyCell = 0.0;
	  _macroFields.Entropy.addArray(cells,EntropyCell);
	  
	  shared_ptr<TArray> EntropyGenRateCell(new TArray(cells.getCount()));
	  *EntropyGenRateCell = 0.0;
	  _macroFields.EntropyGenRate.addArray(cells,EntropyGenRateCell);

	  shared_ptr<TArray> EntropyGenRateColl(new TArray(cells.getCount()));
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
        shared_ptr<TArray> KnqCell(new TArray(cells.getCount()));
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
            if((_bcMap[fg.id]->bcType == "WallBC")||(_bcMap[fg.id]->bcType == "RealWallBC")||(_bcMap[fg.id]->bcType == "SymmetryBC"))
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
            else
	    {
                const StorageSite& faces = fg.site;
                const int faceCount=faces.getCount();
                const int offSet=faces.getOffset();
                         
                for(int i=offSet;i<offSet+faceCount;i++)
                  BCfArray[i]=0;
	    } 
	}
	
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

	      int newCount= MakeCoarseMesh(finerModel->getMeshList(),
					   finerModel->getGeomFields(),
					   *newMeshesPtr);
                    
	      TCOMET* newModelPtr=new COMETModel(*newMeshesPtr,thisLevel,
						 finerModel->getGeomFields(),
						 *newMacroPtr,*newQuadPtr);
              cout<<"Number of cells in level "<<thisLevel<<"  is "<<newCount<<endl;            
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
    
	      newModelPtr->MakeCoarseModel(newModelPtr);
	  }
      }
      else if(_options.AgglomerationMethod=="AMG")
	throw CException("Have not implemented AMG agglomeration method.");
      else
	throw CException("Unknown agglomeration method.");
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

  int MakeCoarseMesh(const MeshList& inMeshes, GeomFields& inGeomFields,
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
        IntArray FineToCoarse(inCellTotal);
        FineToCoarse=-1;

        const CRConnectivity& inCellinFaces=mesh.getCellFaces();
        const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
        Field& areaMagField=inGeomFields.areaMag;
        const TArray& areaMagArray=dynamic_cast<const TArray&>(areaMagField[inFaces]);
        const BCfaceArray& inBCfArray=*(_BFaces[n]);

        //first sweep to make initial pairing
        int pairWith;
        for(int c=0;c<inCellCount;c++)
          {
            if(FineToCoarse[c]<0) //dont bother if im already paired
              {
                //loop through all neighbors to find pairing
                const int neibCount=inCellinFaces.getCount(c);
                pairWith=-1;
                T maxArea=0.;
                int c2;
                for(int neib=0;neib<neibCount;neib++)
                  {
                    const int f=inCellinFaces(c,neib);
                    
                    if(inBCfArray[f]==0)  //not a boundary face
                      {
                        if(c==inFaceinCells(f,1))
                          c2=inFaceinCells(f,0);
                        else
                          c2=inFaceinCells(f,1);

                        if(FineToCoarse[c2]==-1)
                          if(areaMagArray[c2]>maxArea)
                            pairWith=c2;
                      }
                  }

                if(pairWith!=-1)
                  {
                    FineToCoarse[c]=coarseCount;
                    FineToCoarse[c2]=coarseCount;
                    coarseCount++;
                  }
              }
          }

        //second sweep to group stragglers, or group with self
        for(int c=0;c<inCellCount;c++)
          {
            if(FineToCoarse[c]==-1)
              {
                const int neibCount=inCellinFaces.getCount(c);
                T maxArea=0.;
                int c2,c2perm;
                pairWith=-1;

                for(int neib=0;neib<neibCount;neib++)
                  {
                    const int f=inCellinFaces(c,neib);
                    
                    if(inBCfArray[f]==0)  //not a boundary face
                      {
                        if(c==inFaceinCells(f,1))
                          c2=inFaceinCells(f,0);
                        else
                          c2=inFaceinCells(f,1);

                        if(areaMagArray[c2]>maxArea)
                          {
                            pairWith=FineToCoarse[c2]; //coarse level cell
                            c2perm=c2;                 //fine level cell
                          }
                      }
                  }

                if(pairWith==-1)
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

        int coarseGhost=coarseCount;
        for(int c=inCellCount;c<inCellTotal;c++)
          {
            FineToCoarse[c]=coarseGhost;
            coarseGhost++;
          }

        //make the coarse cell to fine cell connectivity.
        outCells.setCount(coarseCount,inGhost);
        CRPtr CoarseToFineCells=CRPtr(new CRConnectivity(outCells,inCells));
        CoarseToFineCells->initCount();

        for(int c=0;c<inCellTotal;c++)
          CoarseToFineCells->addCount(FineToCoarse[c],1);

        CoarseToFineCells->finishCount();

        for(int c=0;c<inCellTotal;c++)
          CoarseToFineCells->add(FineToCoarse[c],c);

        CoarseToFineCells->finishAdd();

        //connectivity between itself (cells) and its finer mesh cells. 
        newMeshPtr->setConnectivity(outCells,inCells,CoarseToFineCells);

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

        outFaces.setCount(counter);

        CRPtr CoarseCellCoarseFace=CRPtr(new CRConnectivity(outCells,outFaces));
        CoarseCellCoarseFace->initCount();

        for(int c=0;c<outCells.getCount();c++)
          {
            const int neibs=CellCellCoarse->getCount(c);
            CoarseCellCoarseFace->addCount(c,neibs);
          }

        CoarseCellCoarseFace->finishCount();

        //make cell connectivity to interior faces.
        IntArray neibCounter(outCells.getCount());
        neibCounter=0;
        counter=0;
        counted=false;
        for(int c=0;c<outCells.getSelfCount();c++)
          {
            counted[c]=true;
            const int neibs=CellCellCoarse->getCount(c);
            for(int n=0;n<neibs;n++)
              {
                const int c1=(*CellCellCoarse)(c,n);
                if(!counted[c1] && c1<outCells.getSelfCount())
                  {
                    CoarseCellCoarseFace->add(c,counter);
                    CoarseCellCoarseFace->add(c1,counter);
                    counter++;
                    neibCounter[c]++;
                    neibCounter[c1]++;
                  }
              }
          }

        //make cell connectivity to boundary faces.
        for(int c=outCells.getSelfCount();c<outCells.getCount();c++)
          {
            const int c1=(*CellCellCoarse)(c,0);
            CoarseCellCoarseFace->add(c1,counter);
            CoarseCellCoarseFace->add(c,counter);
            counter++;
            neibCounter[c]++;
            neibCounter[c1]++;
          }

        CoarseCellCoarseFace->finishAdd();

        CRPtr CoarseFaceCoarseCell=CoarseCellCoarseFace->getTranspose();

        newMeshPtr->setConnectivity(outCells,outFaces,CoarseCellCoarseFace);
        newMeshPtr->setConnectivity(outFaces,outCells,CoarseFaceCoarseCell);

        CRPtr CoarseFacesFineFaces=CRPtr(new CRConnectivity(outFaces,inFaces));
        CoarseFacesFineFaces->initCount();

        for(int f=0;f<inFaceCount;f++)
          {
            int fc0=inFaceinCells(f,0);
            int fc1=inFaceinCells(f,1);
            const int cc0=FineToCoarse[fc0];
            const int cc1=FineToCoarse[fc1];

            if(cc1!=cc0)
              {
                const int cfaces=CoarseCellCoarseFace->getCount(cc0);

                for(int cf=0;cf<cfaces;cf++)
                  {
                    const int face=(*CoarseCellCoarseFace)(cc0,cf);
                    const int tempc0=(*CoarseFaceCoarseCell)(face,0);
                    const int tempc1=(*CoarseFaceCoarseCell)(face,1);
                    
                    if(((cc0==tempc0)&&(cc1==tempc1))||((cc1==tempc0)&&(cc0==tempc1)))
                      {
                        CoarseFacesFineFaces->addCount(face,1);
                        break;
                      }
                  }
              }
          }

        CoarseFacesFineFaces->finishCount();

        for(int f=0;f<inFaceCount;f++)
          {
            int fc0=inFaceinCells(f,0);
            int fc1=inFaceinCells(f,1);
            const int cc0=FineToCoarse[fc0];
            const int cc1=FineToCoarse[fc1];
            if(cc1!=cc0)
              {
                const int cfaces=CoarseCellCoarseFace->getCount(cc0);

                for(int cf=0;cf<cfaces;cf++)
                  {
                    const int face=(*CoarseCellCoarseFace)(cc0,cf);
                    const int tempc0=(*CoarseFaceCoarseCell)(face,0);
                    const int tempc1=(*CoarseFaceCoarseCell)(face,1);
                    
                    if(((cc0==tempc0)&&(cc1==tempc1))||((cc1==tempc0)&&(cc0==tempc1)))
                      {
                        CoarseFacesFineFaces->add(face,f);
                        break;
                      }
                  }
              }
          }

        CoarseFacesFineFaces->finishAdd();

        const int interiorCount=outFaces.getCount()-inGhost;

        const StorageSite& interiorFaces=newMeshPtr->createInteriorFaceGroup(interiorCount);

        int inOffset=interiorCount;
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
          {
            const FaceGroup& fg=*fgPtr;
            const int size=fg.site.getCount();
            newMeshPtr->createBoundaryFaceGroup(size,inOffset,fg.id,fg.groupType);
            inOffset+=size;
          }

        //now make the geom fields
        const int outCellsCount=outCells.getSelfCount();
        TArrptr outCellVolumePtr=TArrptr(new TArray(outCellsCount));
        TArray& outCV=*outCellVolumePtr;
        outCV=0.;

        Field& VolumeField=inGeomFields.volume;
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

        VolumeField.addArray(outCells,outCellVolumePtr);

        const int outFacesCount=outFaces.getCount();
        VT3Ptr outFaceAreaPtr=VT3Ptr(new VectorT3Array(outFacesCount));
        VectorT3Array& outFA=*outFaceAreaPtr;
        TArrptr outFaceAreaMagPtr=TArrptr(new TArray(outFacesCount));
        TArray& outFAMag=*outFaceAreaMagPtr;

        Field& FaceAreaField=inGeomFields.area;
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

        FaceAreaField.addArray(outFaces,outFaceAreaPtr);
        areaMagField.addArray(outFaces,outFaceAreaMagPtr);

	Field& ibTypeField=inGeomFields.ibType;
	shared_ptr<IntArray> ibTypePtr(new IntArray(outCells.getSelfCount()));
	*ibTypePtr = Mesh::IBTYPE_FLUID;
	ibTypeField.addArray(outCells,ibTypePtr);

        if(smallestMesh<0)
          smallestMesh=outCells.getSelfCount();
        else
          {
            if(outCells.getSelfCount()<smallestMesh)
              smallestMesh=outCells.getSelfCount();
          }

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

    return smallestMesh;
  }

  void doSweeps(const int sweeps, const int num)
  {
    for(int sweepNo=0;sweepNo<sweeps;sweepNo++)
      smooth(num);
  }

  void smooth(const int num)
  {
    const int numMeshes=_meshes.size();
    for(int msh=0;msh<numMeshes;msh++)
      {
        const Mesh& mesh=*_meshes[msh];
        const BCcellArray& BCArray=*(_BCells[msh]);
        const BCfaceArray& BCfArray=*(_BFaces[msh]);
	const BCcellArray& ZCArray=*(_ZCells[msh]);
	COMETESBGKDiscretizer<T> CDisc(mesh,_geomFields,_macroFields,_quadrature,
				       _dsfPtr,_dsfPtr1,_dsfPtr2,_dsfEqPtrES,_dsfPtrRes,_dsfPtrFAS,
				       _options["timeStep"],_options.timeDiscretizationOrder,
				       _options.transient,_options["rho_init"], 
				       _options["T_init"],_options["molecularWeight"],
				       _bcMap,_faceReflectionArrayMap,BCArray,BCfArray,ZCArray);

        CDisc.setfgFinder();
        
	CDisc.COMETSolve(1,_level); //forward
	//callCOMETBoundaryConditions();
	ComputeCollisionfrequency();
	//update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
	if (_options.fgamma==0){initializeMaxwellianEq();}
	else{ EquilibriumDistributionBGK();}
	if (_options.fgamma==2){EquilibriumDistributionESBGK();} 
        
	CDisc.COMETSolve(-1,_level); //reverse
	if((num==1)||(num==0&&_level==0))
	{
	    //callCOMETBoundaryConditions();
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
	COMETESBGKDiscretizer<T> CDisc(mesh,_geomFields,_macroFields,_quadrature,
				       _dsfPtr,_dsfPtr1,_dsfPtr2,_dsfEqPtrES,_dsfPtrRes,_dsfPtrFAS,
				       _options["timeStep"],_options.timeDiscretizationOrder,
				       _options.transient,_options["rho_init"], 
				       _options["T_init"],_options["molecularWeight"],
				       _bcMap,_faceReflectionArrayMap,BCArray,BCfArray,ZCArray);

        CDisc.setfgFinder();
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
    doSweeps(_options.preSweeps,1);
    
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
    
    doSweeps(_options.postSweeps,0);
  }

  void injectResid()
  {
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& finerMesh=*_meshes[n];
        const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
        MacroFields& coarserMacro=_coarserLevel->getMacro();
        DistFunctFields<T>& coarserdsf = _coarserLevel->getdsf();
        DistFunctFields<T>& coarserdsf0 = _coarserLevel->getdsf0();
	DistFunctFields<T>& coarserdsfFAS = _coarserLevel->getdsfFAS();
        const StorageSite& finerCells=finerMesh.getCells();
        const StorageSite& coarserCells=coarserMesh.getCells();
        const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);
        const TArray& coarserVol=dynamic_cast<TArray&>(_geomFields.volume[coarserCells]);
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
    cout<<"Initial Residual:"<<_initialResidual<<"  ResidualRatio: "<<residualRatio<<endl;
    int niters=0;
    const T absTol=_options.absoluteTolerance;
    const T relTol=_options.relativeTolerance;
    const int show=_options.showResidual;

    while((niters<iters) && ((_residual>absTol)&&(residualRatio>relTol)))
      {
        cycle();
        niters++;
        _residual=updateResid(false);
	residualRatio=_residual/_initialResidual;
        if(niters%show==0)
          cout<<"Iteration:"<<niters<<" Residual:"<<_residual<<"  ResidualRatio: "<<residualRatio<<endl;

      }
    callCOMETBoundaryConditions();
    //cout<<endl<<"Total Iterations:"<<niters<<" Residual:"<<_residual<<endl;
  }
  
  void  computeSurfaceForce(const StorageSite& solidFaces, bool perUnitArea)
  {
    
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
  Quadrature<T>& _quadrature;
 
  MacroFields& _macroFields;
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
};

#endif
