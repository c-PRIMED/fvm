// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KINETICMODEL_H_
#define _KINETICMODEL_H_

#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

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

#include "KineticBC.h"
#include "TimeDerivativeDiscretization_Kmodel.h"
#include "CollisionTermDiscretization.h"
#include "ConvectionDiscretization_Kmodel.h"
#include "KineticBoundaryConditions.h"
#include "GenericKineticBCS.h"
#include "GenericKineticIBDiscretization.h"

#include "Linearizer.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "CRMatrixTranspose.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"

#include "MatrixOperation.h"
#include "NumType.h"

#include "StressTensor.h"

template<class T>
class KineticModel : public Model
{
 public:
  typedef typename NumTypeTraits<T>:: T_Scalar T_Scalar;
  typedef Array<int> IntArray;
  typedef Array<T> TArray;
  typedef  Array2D<T> TArray2D;
  typedef Vector<T,3> VectorT3; 
  typedef Array<VectorT3> VectorT3Array;
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
  //typedef SquareMatrix<T,5> SqMatrix5;

  typedef std::map<int,KineticBC<T>*> KineticBCMap;
  typedef std::map<int,KineticVC<T>*> KineticVCMap;
  /**
   * Calculation of macro-parameters density, temperature, components of velocity, pressure
   * by taking moments of distribution function using quadrature points and weights from quadrature.h
   */
  //MacroFields& macroFields;
  
  KineticModel(const MeshList& meshes, const GeomFields& geomFields, MacroFields& macroFields, const Quadrature<T>& quad):
       
    Model(meshes),
    _geomFields(geomFields),
    _quadrature(quad),
    _macroFields(macroFields),
    _dsfPtr(_meshes,_quadrature,"dsf_"),
    _dsfPtr1(_meshes,_quadrature,"dsf1_"),
    _dsfPtr2(_meshes,_quadrature,"dsf2_"),
    _dsfEqPtr(_meshes,_quadrature,"dsfEq_"),
    _dsfEqPtrES(_meshes,_quadrature,"dsfEqES_"),
    _initialKmodelNorm(),
    _niters(0)
    {     
     
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  
	  KineticVC<T> *vc(new KineticVC<T>());
	  vc->vcType = "flow";
	  _vcMap[mesh.getID()] = vc;
	}
      init(); 
      SetBoundaryConditions();
      
      //weightedMaxwellian(1.0,0.00,0.00,1.0,1.0);
      InitializeMacroparameters();
      initializeMaxwellian();
      
      ComputeMacroparameters(); //calculate density,velocity,temperature
    
      ComputeCollisionfrequency(); //calculate viscosity, collisionFrequency
      initializeMaxwellianEq();    //equilibrium distribution
      
      //EquilibriumDistributionBGK();
      //callBoundaryConditions();
     
    }
    
    void init()

    {
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  
	  const KineticVC<T>& vc = *_vcMap[mesh.getID()];
	  
	  const StorageSite& cells = mesh.getCells();
	  
	  const int nCells = cells.getCountLevel1();
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
	  *rhoCell = vc["density"];
	  _macroFields.density.addArray(cells,rhoCell);
	  
	  shared_ptr<TArray> muCell(new TArray(nCells));
	  *muCell = vc["viscosity"];
	  _macroFields.viscosity.addArray(cells,muCell);
	  
	  shared_ptr<TArray> tempCell(new TArray(cells.getCountLevel1()));
	  *tempCell = _options["operatingTemperature"];
	  _macroFields.temperature.addArray(cells,tempCell);
	  
	  shared_ptr<TArray> collFreqCell(new TArray(cells.getCountLevel1()));
	  *collFreqCell = vc["viscosity"];
	  _macroFields.collisionFrequency.addArray(cells,collFreqCell);
	  
	  
	  //coeffs for perturbed BGK distribution function
	  shared_ptr<VectorT5Array> coeffCell(new VectorT5Array(cells.getCountLevel1()));
	  VectorT5 initialCoeff;
	  initialCoeff[0] = 1.0;
	  initialCoeff[1] = 1.0;
	  initialCoeff[2] = 0.0; 
	  initialCoeff[3] = 0.0;
	  initialCoeff[4] = 0.0;
	  *coeffCell = initialCoeff;
	  _macroFields.coeff.addArray(cells,coeffCell);
	  
	  //coeffs for perturbed BGK distribution function
	  shared_ptr<VectorT10Array> coeffgCell(new VectorT10Array(cells.getCountLevel1()));
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

	const int numDirections = _quadrature.getDirCount();
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	//FILE * pFile;
	//pFile=fopen("ref_incMEMOSA.txt","w");
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups()){
	  const FaceGroup& fg = *fgPtr; 
	  
	  if((fg.groupType == "symmetry")||(fg.groupType == "realwall")||(fg.groupType == "velocity-inlet")){
	  
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
	_niters  =0;
	_initialKmodelNorm = MFRPtr();
	//_initialKmodelvNorm = MFRPtr();
	
      }
  }
    
  void InitializeMacroparameters()
  {  const int numMeshes =_meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1(); 
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
	const int nCells = cells.getCountLevel1(); 
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
	const int nCells = cells.getCountLevel1();  //
	
	
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

  
  
  void ComputeMacroparametersESBGK() 
  {  
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {	
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();
	
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
	const int nCells = cells.getCountLevel1();
	
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
	const int nCells = cells.getCountLevel1();
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
	const int nCells = cells.getCountLevel1();
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
	const int nCells = cells.getCountLevel1();
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
	const int nCells = cells.getCountLevel1();
	
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
	const int nCells = cells.getCountLevel1();
	
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
	const int nCells = cells.getCountLevel1();
	
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
	const int nCells = cells.getCountLevel1();


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
	const int nCells = cells.getCountLevel1();

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

  void weightedMaxwellian(double weight1,double uvel1,double vvel1,double wvel1,double uvel2,double vvel2,double wvel2,double temp1,double temp2)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCountLevel1();
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
	const int nCells = cells.getCountLevel1();
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
  KineticBCMap& getBCMap()  {return _bcMap;}
  KineticVCMap& getVCMap()  {return _vcMap;}
  
  KineticModelOptions<T>&   getOptions() {return _options;}
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
	
	KineticVC<T> *vc(new KineticVC<T>());
	vc->vcType = "flow";
	_vcMap[mesh.getID()] = vc;
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if (_bcMap.find(fg.id) == _bcMap.end())
	      {
		KineticBC<T> *bc(new KineticBC<T>());
		
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
		  throw CException("KineticModel: unknown face group type "
				     + fg.groupType);
	      }
	  }
	/*
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if (_bcMap.find(fg.id) == _bcMap.end())
	      { KineticBC<T> *bc(new KineticBC<T>());
		
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
  
  
  void initKineticModelLinearization(LinearSystem& ls, int direction)
  {
   const int numMeshes = _meshes.size();
   for (int n=0; n<numMeshes; n++)
     {
       const Mesh& mesh = *_meshes[n];
       const StorageSite& cells = mesh.getCells();
       
       Field& fnd = *_dsfPtr.dsf[direction];
       //const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
       
       MultiField::ArrayIndex vIndex(&fnd,&cells); //dsf is in 1 direction
       
       ls.getX().addArray(vIndex,fnd.getArrayPtr(cells));
       
       const CRConnectivity& cellCells = mesh.getCellCells();
       
       shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells)); //diagonal off diagonal, phi DiagTensorT3,T,VectorT3 for mometum in fvmbase 
       
       ls.getMatrix().addMatrix(vIndex,vIndex,m);
     }
   //mesh.getBoudaryfaces and interfacefaces?
  } 
  
  
  void linearizeKineticModel(LinearSystem& ls, int direction)
  {
    // _velocityGradientModel.compute();
    DiscrList discretizations;
    
   
    
    Field& fnd = *_dsfPtr.dsf[direction]; 
    //Field& feq = *_dsfEqPtr.dsf[direction];
    //if(_options.ESBGK_fgamma){feq = *_dsfEqPtrES.dsf[direction];}
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
   
    
    if(_options.fgamma==2){ 
      Field& feqES = *_dsfEqPtrES.dsf[direction];
      shared_ptr<Discretization>
	sdEQ(new CollisionTermDiscretization<T,T,T>
	     (_meshes, _geomFields, 
	      fnd,feqES,
	      _macroFields.collisionFrequency)); 
      discretizations.push_back(sdEQ);}
    else{
      Field& feq = *_dsfEqPtr.dsf[direction];
      shared_ptr<Discretization>
	sd(new CollisionTermDiscretization<T,T,T>
	   (_meshes, _geomFields, 
	    fnd,feq,
	    _macroFields.collisionFrequency));
    discretizations.push_back(sd);
    } 
    
    
    shared_ptr<Discretization> 
      cd(new ConvectionDiscretization_Kmodel<T,T,T> 
	 (_meshes,
	  _geomFields,
	  fnd,
	  cx[direction],
	  cy[direction],
 	  cz[direction],
	  _options.CentralDifference
	  //_options["nonDimLt"],
	  //_options["nonDimLx"],_options["nonDimLy"],_options["nonDimLz"],
	  ));
    discretizations.push_back(cd);
    
    if (_options.transient)
      {
	// const int direction(0);  
	Field& fnd = *_dsfPtr.dsf[direction];            
	Field& fnd1 = *_dsfPtr1.dsf[direction];
	Field& fnd2 = *_dsfPtr2.dsf[direction];

        
	shared_ptr<Discretization>
	  td(new TimeDerivativeDiscretization_Kmodel<T,T,T>
	     (_meshes,_geomFields,
	      fnd,fnd1,fnd2,
	      _options["timeStep"],
	      _options["nonDimLt"],
	      _options.timeDiscretizationOrder));
	
	discretizations.push_back(td);
	
      }
    if(_options.ibm_enable){ 
    shared_ptr<Discretization>
      ibm(new GenericKineticIBDiscretization<T,T,T>
	  (_meshes,_geomFields,fnd,
	   cx[direction],
	   cy[direction],
	   cz[direction],
	   _macroFields));
    discretizations.push_back(ibm);
    }

    Linearizer linearizer;
    
    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
			 ls.getX(), ls.getB());
    
    // boundary conditions
    const double epsilon=_options.epsilon_ES;
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();

	    const KineticBC<T>& bc = *_bcMap[fg.id];
	    Field& fnd = *_dsfPtr.dsf[direction]; //field in a direction
	    TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
	    BaseGenericKineticBCS<T,T,T> gkbc(faces, mesh, _geomFields,
						 fnd,
						 ls.getMatrix(),
						 ls.getX(),
					      ls.getB());
					      
	   
	    const CRConnectivity& faceCells = mesh.getFaceCells(faces);			 
	    const Field& areaMagField = _geomFields.areaMag;
            const TArray& faceAreaMag = dynamic_cast<const TArray &>(areaMagField[faces]);
	    const Field& areaField = _geomFields.area;
	    const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]); 
	     
	    FloatValEvaluator<VectorT3>
	      bVelocity(bc.getVal("specifiedXVelocity"),
			bc.getVal("specifiedYVelocity"),
			bc.getVal("specifiedZVelocity"),
			faces);

	    if (( bc.bcType == "ZeroGradBC")) 
	      {
		for(int f=0; f< nFaces; f++)
		  {const int c1= faceCells(f,1);// boundary cell
		     T bvalue =dsf[c1];
		      gkbc.applyDirichletBC(f,bvalue);
		      // gkbc.applyExtrapolationBC(f);
		  }
	      } 
	    else if (( bc.bcType == "VelocityInletBC")){
	      for(int f=0; f< nFaces; f++)
		{
		  const VectorT3 en = faceArea[f]/faceAreaMag[f];
		  const T c_dot_en = cx[direction]*en[0]+cy[direction]*en[1]+cz[direction]*en[2];		 
		  if(c_dot_en  < T_Scalar(epsilon))
		    //incoming direction - dirchlet bc
		    { const int c1= faceCells(f,1);
		      T bvalue =dsf[c1];
		      gkbc.applyDirichletBC(f,bvalue);
		    }
		  else{
		    //outgoing direction - extrapolation bc
		    gkbc.applyExtrapolationBC(f);
		  } 
		}
	    }
	    //if ((bc.bcType == "WallBC")||(bc.bcType=="PressureInletBC")|| (bc.bcType=="PressureOutletBC")|| (bc.bcType=="VelocityInletBC"))
	   
	    else{
	      for(int f=0; f< nFaces; f++)
		{
		  const VectorT3 en = faceArea[f]/faceAreaMag[f];
		  const T c_dot_en = cx[direction]*en[0]+cy[direction]*en[1]+cz[direction]*en[2];
		  const VectorT3  WallVelocity = bVelocity[f];
		  const T uwall = WallVelocity[0]; const T vwall = WallVelocity[1];
		  const T wwall = WallVelocity[2]; const T wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
		  if(c_dot_en -wallV_dot_en < T_Scalar(epsilon))
		    //incoming direction - dirchlet bc
		    { const int c1= faceCells(f,1);// boundary cell
		      T bvalue =dsf[c1];
		      gkbc.applyDirichletBC(f,bvalue);
		    }
		  else{
		    //outgoing direction - extrapolation bc
		    gkbc.applyExtrapolationBC(f);
		  } 
		  
		}
	    }
  
	    
	    
	  }

	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr; 
	    const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
	    //const KineticBC<T>& bc = *_bcMap[fg.id];
	    //Field& fnd = *_dsfPtr.dsf[direction]; //field in a direction
	    //TArray& dsf = dynamic_cast< TArray&>(fnd[cells]);
	    BaseGenericKineticBCS<T,T,T> gkbc(faces, mesh, _geomFields,
						 fnd,
						 ls.getMatrix(),
						 ls.getX(),
					      ls.getB());
	    for(int f=0; f< nFaces; f++)
		 {gkbc.applyInterfaceBC(f);}//do nothign
	    
	  }


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
	   

           shared_ptr<TArray> ibV(new TArray(ibFaces.getCount()));
       
           const TArray& cV =
            dynamic_cast<const TArray&>(fnd[cells]);

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

          fnd.addArray(ibFaces,ibV);	 
	    }
	  }
      }
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
	if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& ibFaces = mesh.getIBFaces();
        
	  GeomFields::SSPair key1(&ibFaces,&cells);
	  const IMatrix& mIC =
	    dynamic_cast<const IMatrix&>
	    (*_geomFields._interpolationMatrices[key1]);
	      
	  IMatrix mICV(mIC);
	  IMatrixV3 mICV3(mIC);  

	  GeomFields::SSPair key2(&ibFaces,&solidFaces);
	  const IMatrix& mIP =
	    dynamic_cast<const IMatrix&>
	    (*_geomFields._interpolationMatrices[key2]);

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

	   //nue interpolation (cells)
	  mICV.multiplyAndAdd(*ibVnue,cNue);
	  mIPV.multiplyAndAdd(*ibVnue,sNue);
	  _macroFields.collisionFrequency.addArray(ibFaces,ibVnue);
	  //temperature interpolation (cells+solidfaces)         
	  mICV.multiplyAndAdd(*ibVtemp,cTemp);
	  mIPV.multiplyAndAdd(*ibVtemp,sTemp);
	  _macroFields.temperature.addArray(ibFaces,ibVtemp);
	  //density interpolation (cells+solidfaces)         
	  mICV.multiplyAndAdd(*ibVdensity,cDensity);
	  mIPV.multiplyAndAdd(*ibVdensity,sDensity);
	  _macroFields.density.addArray(ibFaces,ibVdensity);
	  //velocity interpolation (cells+solidfaces) 
	  mICV3.multiplyAndAdd(*ibVvel,cVel);
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

//#ifdef FVM_PARALLEL
//	      MPI::COMM_WORLD.Allreduce( MPI::IN_PLACE,dsf.getData(),solidFaces.getCount() , MPI::DOUBLE, MPI::SUM);
//#endif
	      const Mesh& mesh = *_meshes[n];
	      if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
		const StorageSite& ibFaces = mesh.getIBFaces();
		TArray& dsfIB = dynamic_cast< TArray&>(fnd[ibFaces]);
		Field& fndEqES = *_dsfEqPtrES.dsf[j];
		TArray& dsfEqES = dynamic_cast< TArray&>(fndEqES[ibFaces]);
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

		for(int f=0; f<nibFaces; f++)
		  {
		    dsfIB[f]=0.0;				  
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
	    if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	    
	      const StorageSite& cells = mesh.getCells();
	      
	    
	      GeomFields::SSPair key1(&solidFaces,&cells);
	      const IMatrix& mIC =
		dynamic_cast<const IMatrix&>
		(*_geomFields._interpolationMatrices[key1]);
	      
	      IMatrix mICV(mIC);
	      const TArray& cV =
		dynamic_cast<const TArray&>(fnd[cells]);

	      ibV->zero();        
	     
	      mICV.multiplyAndAdd(*ibV,cV);
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
          
	const int numFields= _quadrature.getDirCount();
	for (int direction = 0; direction < numFields; direction++)
	  {
	    Field& fnd = *_dsfPtr.dsf[direction];
	    Field& fndEqES = *_dsfEqPtrES.dsf[direction];
	    const Mesh& mesh = *_meshes[n];
	    if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
	
	      const StorageSite& cells = mesh.getCells();
	      const StorageSite& faces = mesh.getFaces();
	      const StorageSite& ibFaces = mesh.getIBFaces();
	
	      GeomFields::SSPair key1(&faces,&cells);
	      const IMatrix& mIC =
		dynamic_cast<const IMatrix&>
		(*_geomFields._interpolationMatrices[key1]);
	      
	      IMatrix mICV(mIC);
	  
	      const TArray& cf =
		dynamic_cast<const TArray&>(fnd[cells]);
	      const TArray& cfEq =
		dynamic_cast<const TArray&>(fndEqES[cells]);

	      shared_ptr<TArray> ibVf(new TArray(ibFaces.getCount()));
	      ibVf->zero();

	      if (_options.fgamma==2){
		shared_ptr<TArray> ibVfEq(new TArray(ibFaces.getCount()));
		ibVfEq->zero();
		mICV.multiplyAndAdd(*ibVfEq,cfEq);
		fndEqES.addArray(ibFaces,ibVfEq);
	      }
	      mICV.multiplyAndAdd(*ibVf,cf);
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
	  if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

	    const StorageSite& cells = mesh.getCells();
	    const StorageSite& ibFaces = mesh.getIBFaces();
        
	  GeomFields::SSPair key1(&ibFaces,&cells);
	  const IMatrix& mIC =
	    dynamic_cast<const IMatrix&>
	    (*_geomFields._interpolationMatrices[key1]);
	      
	  IMatrix mICV(mIC);
	  IMatrixV3 mICV3(mIC);  

	  GeomFields::SSPair key2(&ibFaces,&solidFaces);
	  const IMatrix& mIP =
	    dynamic_cast<const IMatrix&>
	    (*_geomFields._interpolationMatrices[key2]);

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

	   //nue interpolation (cells)
           mICV.multiplyAndAdd(*ibVnue,cNue);
	   mIPV.multiplyAndAdd(*ibVnue,sNue);
           _macroFields.collisionFrequency.addArray(ibFaces,ibVnue);
	   //temperature interpolation (cells+solidfaces)         
	   mICV.multiplyAndAdd(*ibVtemp,cTemp);
	   mIPV.multiplyAndAdd(*ibVtemp,sTemp);
           _macroFields.temperature.addArray(ibFaces,ibVtemp);
	   //density interpolation (cells+solidfaces)         
	   mICV.multiplyAndAdd(*ibVdensity,cDensity);
	   mIPV.multiplyAndAdd(*ibVdensity,sDensity);
           _macroFields.density.addArray(ibFaces,ibVdensity);
	   //velocity interpolation (cells+solidfaces) 
	   mICV3.multiplyAndAdd(*ibVvel,cVel);
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
		if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){
		  const StorageSite& ibFaces = mesh.getIBFaces();
		  const CRConnectivity& solidFacesToibFaces
		    = mesh.getConnectivity(solidFaces,ibFaces);
		  const IntArray& ibFaceIndices = mesh.getIBFaceList();
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
		      const StorageSite& faces = mesh.getFaces();
		      const int c = sFCCol[nc];
		      const int faceIB= ibFaceIndices[c];
		      const VectorT3Array& solidFaceCentroid =
			dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[solidFaces]);
		      const VectorT3Array& faceCentroid =
			dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[faces]);
	      
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
		      const StorageSite& faces = mesh.getFaces();
		      const VectorT3Array& solidFaceCentroid =
			dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[solidFaces]);
		      const VectorT3Array& faceCentroid =
			dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[faces]);
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

 void ConservationofMFSolid(const StorageSite& solidFaces, const int output =0, bool perUnitArea=0) const
 {
   FILE * pFile;
	    
   shared_ptr<VectorT6Array> Stress(new VectorT6Array(solidFaces.getCount()));
   Stress->zero();
   _macroFields.Stress.addArray(solidFaces,Stress);

   shared_ptr<VectorT3Array> Force(new VectorT3Array(solidFaces.getCount()));
   Force->zero();
   _macroFields.force.addArray(solidFaces,Force);
	        
   shared_ptr<TArray> TemperatureIB(new TArray(solidFaces.getCount()));
   TemperatureIB->zero();
   _macroFields.temperatureIB.addArray(solidFaces,TemperatureIB);

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
	  dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[solidFaces]);
	const VectorT3Array& solidFaceArea =
	  dynamic_cast<const VectorT3Array&>(_geomFields.area[solidFaces]);
	const TArray& solidFaceAreaMag =
	  dynamic_cast<const TArray&>(_geomFields.areaMag[solidFaces]);
	const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[solidFaces]);
	TArray& density  = dynamic_cast<TArray&>(_macroFields.density[solidFaces]);
	TArray& temperature  = dynamic_cast<TArray&>(_macroFields.temperature[solidFaces]);
	TArray& tempIB  = dynamic_cast<TArray&>(_macroFields.temperatureIB[solidFaces]);
	VectorT3Array& force  = dynamic_cast<VectorT3Array&>(_macroFields.force[solidFaces]);
	VectorT6Array& stress  = dynamic_cast<VectorT6Array&>(_macroFields.Stress[solidFaces]);
	const T Lx=_options["nonDimLx"];
	const T Ly=_options["nonDimLy"];
	const T Lz=_options["nonDimLz"];
    
	const T uwall = v[i][0];
	const T vwall = v[i][1];
	const T wwall = v[i][2];
	const T Twall = temperature[i];

	T Nmr(0.0) ;
	T Dmr(0.0) ;
	T incomFlux(0.0);
	T TempWall(0.0);
	T mWall(0.0);

	for (int j=0; j<numDirections; j++)
	  {		
	    Field& fnd = *_dsfPtr.dsf[j];
	    TArray& dsf = dynamic_cast< TArray&>(fnd[solidFaces]);
	    const VectorT3 en = solidFaceArea[i]/solidFaceAreaMag[i];
	    const T c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	    const T wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	    const T fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	    
	    if (c_dot_en-wallV_dot_en > 0) //outgoing
	      {
		Dmr = Dmr - fwall*wts[j]*(c_dot_en-wallV_dot_en);
		incomFlux=incomFlux-dsf[i]*wts[j]*(c_dot_en-wallV_dot_en);
	      }
	    else
	      {
		Nmr = Nmr + dsf[i]*wts[j]*(c_dot_en-wallV_dot_en);
          	if(output==1){
               TempWall = TempWall + dsf[i]*(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))*wts[j]*0.5;
                density[i] = density[i] + dsf[i]*wts[j]*0.5;
		stress[i][0] +=pow((cx[j]-uwall),2.0)*dsf[i]*wts[j]*0.5;
		stress[i][1] +=pow((cy[j]-vwall),2.0)*dsf[i]*wts[j]*0.5;
	        stress[i][2] +=pow((cz[j]-wwall),2.0)*dsf[i]*wts[j]*0.5;
		stress[i][3] +=(cx[j]-uwall)*(cy[j]-vwall)*dsf[i]*wts[j]*0.5;
		stress[i][4] +=(cy[j]-vwall)*(cz[j]-wwall)*dsf[i]*wts[j]*0.5;
		stress[i][5] +=(cx[j]-uwall)*(cz[j]-wwall)*dsf[i]*wts[j]*0.5;
                }
	      }
	  }
	const T nwall = Nmr/Dmr; // incoming wall number density for initializing Maxwellian
		    	    
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
 	        if(output==1){
                TempWall = TempWall + dsf[i]*(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))*wts[j]*0.5;
                density[i] = density[i] + dsf[i]*wts[j]*0.5;
		stress[i][0] +=pow((cx[j]-uwall),2.0)*dsf[i]*wts[j]*0.5;
		stress[i][1] +=pow((cy[j]-vwall),2.0)*dsf[i]*wts[j]*0.5;
	        stress[i][2] +=pow((cz[j]-wwall),2.0)*dsf[i]*wts[j]*0.5;
		stress[i][3] +=(cx[j]-uwall)*(cy[j]-vwall)*dsf[i]*wts[j]*0.5;
		stress[i][4] +=(cy[j]-vwall)*(cz[j]-wwall)*dsf[i]*wts[j]*0.5;
		stress[i][5] +=(cx[j]-uwall)*(cz[j]-wwall)*dsf[i]*wts[j]*0.5;
               }
	      }	    
	    else
	      dsf[i]=dsf[i];
	  }
	  if(output==1){
	  const VectorT3& Af = solidFaceArea[i];
	  force[i][0] = Af[0]*Ly*Lz*stress[i][0] + Af[1]*Lz*Lx*stress[i][3] + Af[2]*Lx*Ly*stress[i][5];
	  force[i][1] = Af[0]*Ly*Lz*stress[i][3] + Af[1]*Lz*Lx*stress[i][1] + Af[2]*Lx*Ly*stress[i][4];
	  force[i][2] = Af[0]*Ly*Lz*stress[i][5] + Af[1]*Lz*Lx*stress[i][4] + Af[2]*Ly*Ly*stress[i][2];
	 
	  pFile = fopen("WallTemperature.dat","a");	  
	   fprintf(pFile,"%E %E %E %E\n",solidFaceCentroid[i][0],solidFaceCentroid[i][1],solidFaceCentroid[i][2],TempWall);
          fclose(pFile);
	}
      }
 }

void MacroparameterIBCell(const StorageSite& solidFaces) const
 {
//    typedef CRMatrixTranspose<T,T,T> IMatrix;
//    typedef CRMatrixTranspose<T,VectorT3,VectorT3> IMatrixV3;
    const int numMeshes = _meshes.size();
//    for (int n=0; n<numMeshes; n++)
//      {
//	const Mesh& mesh = *_meshes[n];
//	if (!mesh.isShell() && mesh.getIBFaces().getCount() > 0){

//	  const StorageSite& cells = mesh.getCells();
//	  const StorageSite& ibFaces = mesh.getIBFaces();
        
//	  GeomFields::SSPair key1(&ibFaces,&cells);
//	  const IMatrix& mIC =
//	    dynamic_cast<const IMatrix&>
//	    (*_geomFields._interpolationMatrices[key1]);
	      
//	  IMatrix mICV(mIC);
//	  IMatrixV3 mICV3(mIC);  

//	  GeomFields::SSPair key2(&ibFaces,&solidFaces);
//	  const IMatrix& mIP =
//	    dynamic_cast<const IMatrix&>
//	    (*_geomFields._interpolationMatrices[key2]);

//	  IMatrix mIPV(mIP);
//	  IMatrixV3 mIPV3(mIP);    

//	  shared_ptr<TArray> ibVtemp(new TArray(ibFaces.getCount()));
	  
//	  const TArray& cTemp  = 
//	    dynamic_cast<TArray&>(_macroFields.temperature[cells]);
//	  const TArray& sTemp  = 
//	    dynamic_cast<TArray&>(_macroFields.temperatureIB[solidFaces]);
	  
//	  ibVtemp->zero();

	  //temperature interpolation (cells+solidfaces)         
//	  mICV.multiplyAndAdd(*ibVtemp,cTemp);
//	  mIPV.multiplyAndAdd(*ibVtemp,sTemp);
//	  _macroFields.temperature.addArray(ibFaces,ibVtemp);
//	}
//      }
    for (int n=0;n<numMeshes;n++)
     {
    const Mesh& mesh = *_meshes[n];
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const StorageSite& ibFaces = mesh.getIBFaces();
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getCountLevel1();
    const StorageSite& faces = mesh.getFaces();
//    const TArray& TempIB  = 
//       dynamic_cast<TArray&>(_macroFields.temperature[ibFaces]);
    TArray& TempB  = 
       dynamic_cast<TArray&>(_macroFields.pressure[cells]);
    const IntArray& ibType = dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
    const IntArray& ibFaceIndex = dynamic_cast<const IntArray&>(_geomFields.ibFaceIndex[faces]);
 
    TArray xB(nCells);
    TArray wB(nCells);

    xB.zero();
    wB.zero();
      
    const int nFaces = faces.getCount();
    for(int f=0; f<nFaces; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);

        if (((ibType[c0] == Mesh::IBTYPE_FLUID) && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID) && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
            // this is an iBFace, determine which cell is interior and which boundary

            const int ibFace = ibFaceIndex[f];
            if (ibFace < 0)
              throw CException("invalid ib face index");
//            const T xFace = TempIB[ibFace];

            if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
		  xB[c1] = xB[c0];
		  wB[c1]++;
	    }
	
            else
	      {
		  xB[c0] = xB[c1];
		  wB[c0]++;
	      }
	}
        else if ((ibType[c0] == Mesh::IBTYPE_FLUID) &&
            (ibType[c1] == Mesh::IBTYPE_FLUID))
        {
            // leave as  is
        }
        else
        {
            // leave as  is
        }
    }

    // set the phi for boundary cells as average of the ib face values
    for(int c=0; c<nCells; c++)
    {
        if (wB[c] > 0)
	cout << "wb value" << wB[c];
          TempB[c] =  xB[c] / T_Scalar(wB[c]);
	  
    }
 
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
	  

  void callBoundaryConditions()
    
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
	    
	    const KineticBC<T>& bc = *_bcMap[fg.id];
	    
	    KineticBoundaryConditions<T,T,T> kbc(faces, mesh,_geomFields,_quadrature,_macroFields,_dsfPtr);
	    
	    FloatValEvaluator<VectorT3> bVelocity(bc.getVal("specifiedXVelocity"),
						  bc.getVal("specifiedYVelocity"),
						  bc.getVal("specifiedZVelocity"),
						  faces);
	    //	   FloatValEvaluator<StressTensor<T>>
	    FloatValEvaluator<VectorT3> bStress(bc.getVal("specifiedTauxx"),bc.getVal("specifiedTauyy"),
						bc.getVal("specifiedTauzz"),faces);
	    //bc.getVal("specifiedTauxy"),
	    //				      bc.getVal("specifiedTauyz"),bc.getVal("specifiedTauzx"),faces);
	    
	    FloatValEvaluator<T> mdot(bc.getVal("specifiedMassFlowRate"),faces);
	    FloatValEvaluator<T> bTemperature(bc.getVal("specifiedTemperature"),faces);
	    FloatValEvaluator<T> bPressure(bc.getVal("specifiedPressure"),faces);
	    FloatValEvaluator<T> accomCoeff(bc.getVal("accommodationCoefficient"),faces);
	    if (bc.bcType == "WallBC")
	      {			
		kbc.applyDiffuseWallBC(bVelocity,bTemperature);
	      }
	    else if (bc.bcType == "RealWallBC")
	      {
		//kbc.applyRealWallBC(bVelocity,bTemperature,accomCoeff);
		map<int, vector<int> >::iterator pos = _faceReflectionArrayMap.find(fg.id);
		const vector<int>& vecReflection=(*pos).second;
		kbc.applyRealWallBC(bVelocity,bTemperature,accomCoeff,vecReflection);
	      }
	    else if(bc.bcType=="SymmetryBC")
	      {
		//kbc.applySpecularWallBC(); //old boundary works only for cartesian-type quadrature
		map<int, vector<int> >::iterator pos = _faceReflectionArrayMap.find(fg.id);
		const vector<int>& vecReflection=(*pos).second;
		kbc.applySpecularWallBC(vecReflection);
	      } 
	    else if(bc.bcType=="ZeroGradBC")
	      {
		kbc.applyZeroGradientBC();
	      }
	    else if(bc.bcType=="PressureInletBC")
	      {
		kbc.applyPressureInletBC(bTemperature,bPressure);
	      } 
	   
	    else if(bc.bcType=="VelocityInletBC")
	      {	
		map<int, vector<int> >::iterator pos = _faceReflectionArrayMap.find(fg.id);
		const vector<int>& vecReflection=(*pos).second;
		kbc.applyInletBC(bVelocity,bTemperature,mdot,vecReflection);
	      }
	    else if(bc.bcType=="PressureOutletBC")
	      {
	    	kbc.applyPressureOutletBC(bTemperature,bPressure);
	      }
	    
	  }
	foreach(const FaceGroupPtr igPtr, mesh.getInterfaceGroups())
	  {
	    
	    const FaceGroup& ig = *igPtr;
	    const StorageSite& faces = ig.site;
	    //const int nFaces = faces.getCount();
	    
	    //const KineticBC<T>& bc = *_bcMap[fg.id];
	    
	    KineticBoundaryConditions<T,T,T> kbc(faces, mesh,_geomFields,_quadrature,_macroFields,_dsfPtr); 
	    if(ig.groupType=="NSinterface")
	      {
		kbc.applyNSInterfaceBC();//bTemperature,bPressure,bVelocity,bStress);
	      } 
	  }
      }//end of loop through meshes
  }
  
    bool advance(const int niter,const int updated=0)
  {
    for(int n=0; n<niter; n++)  
      {
	const int N123 =_quadrature.getDirCount();
	MFRPtr rNorm;
	//MFRPtr vNorm;
       	//const TArray& cx= dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	//const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

	//callBoundaryConditions();
   	
	/*
	_macroFields.velocity.syncLocal();
	_macroFields.temperature.syncLocal();
	_macroFields.density .syncLocal();
	*/
	for(int direction=0; direction<N123;direction++)
	  {
	    LinearSystem ls;
	    initKineticModelLinearization(ls, direction);
	    ls.initAssembly();
	    linearizeKineticModel(ls,direction);
	    ls.initSolve();
	    
	    //const T* kNorm=_options.getKineticLinearSolver().solve(ls);
	     MFRPtr kNorm(_options.getKineticLinearSolver().solve(ls));
	    
	     
	     if (!rNorm)
	       rNorm = kNorm;
	     else
             {
                 // find the array for the 0the direction residual and
                 // add the current residual to it
                 Field& fn0 = *_dsfPtr.dsf[0];
                 Field& fnd = *_dsfPtr.dsf[direction];
                 
                 ArrayBase& rArray0 = (*rNorm)[fn0];
                 ArrayBase& rArrayd = (*kNorm)[fnd];//*wts[direction];
                 rArray0 += rArrayd;//*wts[direction];
		 
		 // ArrayBase& vArray0 = (*vNorm)[fn0];
		 // vArray0 += rArrayd;//*cx[direction]*wts[direction];
             }

	     ls.postSolve();
	     ls.updateSolution();
	    		 
	     _options.getKineticLinearSolver().cleanup();
	    
	  }

	
	
	if (!_initialKmodelNorm||updated==1) _initialKmodelNorm = rNorm;
	//if (!_initialKmodelvNorm) _initialKmodelvNorm = vNorm;
	if (_niters < 5||updated==2)
        {
             _initialKmodelNorm->setMax(*rNorm);
	     // _initialKmodelvNorm->setMax(*vNorm);
            
        } 
 
	MFRPtr normRatio((*rNorm)/(*_initialKmodelNorm));	
	//	MFRPtr vnormRatio((*vNorm)/(*_initialKmodelvNorm));
#ifdef FVM_PARALLEL
	if ( MPI::COMM_WORLD.Get_rank() == 0 ){
	  if (_options.printNormalizedResiduals){
	    cout << _niters << ": " << *normRatio << endl;}
	  else{
	    cout << _niters << ": " << *rNorm <<endl; }}
#endif
#ifndef FVM_PARALLEL 
	if (_options.printNormalizedResiduals){
	  cout << _niters << ": " << *normRatio << endl;}
	else{
	  cout << _niters << ": " << *rNorm <<endl; }
#endif
	_niters++;
	//break here

	callBoundaryConditions();
	//cout << "called boundary"<<endl;
	ComputeMacroparameters();	//update macroparameters
        ComputeCollisionfrequency();
	
	//update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
	if (_options.fgamma==0){initializeMaxwellianEq();}
	else{ EquilibriumDistributionBGK();}
	//cout << "called BGk" <<endl;
	if (_options.fgamma==2){EquilibriumDistributionESBGK();}


	if ((*rNorm < _options.absoluteTolerance)||(*normRatio < _options.relativeTolerance )){
	  //&& ((*vNorm < _options.absoluteTolerance)||(*vnormRatio < _options.relativeTolerance )))
	  return true;
	}
	
	return false;
	
	
      }
    

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
	    GeomFields::SSPair key1(&solidFaces,&cells);
	    const IMatrix& mIC =
	      dynamic_cast<const IMatrix&>
	      (*_geomFields._interpolationMatrices[key1]);
	    const Array<T>& iCoeffs = mIC.getCoeff();
	    for(int j=0;j<N123;j++){
	      Field& fnd = *_dsfPtr.dsf[j];
	      const TArray& f_dsf = dynamic_cast<const TArray&>(fnd[cells]);
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
 
  const DistFunctFields<T>& getdsf() const { return _dsfPtr;} 
  const DistFunctFields<T>& getdsf1() const { return _dsfPtr1;} 
  const DistFunctFields<T>& getdsf2() const { return _dsfPtr2;}
  const DistFunctFields<T>& getdsfEq() const { return _dsfEqPtr;} 
  const DistFunctFields<T>& getdsfEqES() const { return _dsfEqPtrES;}
    

 private:
  //shared_ptr<Impl> _impl;
 
  const GeomFields& _geomFields;
  const Quadrature<T>& _quadrature;
 
  MacroFields& _macroFields;
  DistFunctFields<T> _dsfPtr;  
  DistFunctFields<T> _dsfPtr1;
  DistFunctFields<T> _dsfPtr2;
  DistFunctFields<T> _dsfEqPtr;
  DistFunctFields<T> _dsfEqPtrES;

  KineticBCMap _bcMap;
  KineticVCMap _vcMap;

  KineticModelOptions<T> _options;
  

  MFRPtr _initialKmodelNorm;
  //MFRPtr _initialKmodelvNorm;
  int _niters;
  map<int, vector<int> > _faceReflectionArrayMap;  
  map<string,shared_ptr<ArrayBase> > _persistenceData;
};

#endif
