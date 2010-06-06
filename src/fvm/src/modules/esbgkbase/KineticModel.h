#ifndef _KINETICMODEL_H_
#define _KINETICMODEL_H_

#include <stdio.h>
#include <math.h>
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

#include "Linearizer.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"



template<class T>
class KineticModel : public Model
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3; 
  typedef Array<VectorT3> VectorT3Array;
  typedef  std::vector<Field*> stdVectorField;
  typedef  DistFunctFields<T> TDistFF;
  
  typedef std::map<int,FlowBC<T>*> FlowBCMap;
  typedef std::map<int,FlowVC<T>*> FlowVCMap;
  /**
   * Calculation of macro-parameters density, temperature, components of velocity, pressure
   * by taking moments of distribution function using quadrature points and weights from quadrature.h
   */
  //MacroFields& macroFields;
  
 KineticModel(const MeshList& meshes, const GeomFields& geomFields, MacroFields& macroFields, const Quadrature<T>& quad):
  //KineticModel(const MeshList& meshes, FlowFields& ffields, const Quadrature<T>& quad):
  
  Model(meshes),
   
    _geomFields(geomFields),
    _quadrature(quad),
    _macroFields(macroFields),
    _dsfPtr(_meshes,_quadrature),
    _dsfPtr1(_meshes,_quadrature),
    _dsfPtr2(_meshes,_quadrature),
    _dsfEqPtr(_meshes,_quadrature)
    {     
      //dsfPtr=new TDistFF(mesh,macroPr,quad);   
      //dsfPtr = new TDistFF(_meshes, quad);     
      // Impl();
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  
	  FlowVC<T> *vc(new FlowVC<T>());
	  vc->vcType = "flow";
	  _vcMap[mesh.getID()] = vc;
	}
      init(); 
      SetBoundaryConditions();
      weightedMaxwellian(0.25,4.0,1.0);

      ComputeMacroparameters(); //calculate density,velocity,temperature
      ComputeCollisionfrequency(); //calculate viscosity, collisionFrequency
      initializeMaxwellianEq();
	// char * filename="f0.plt";
	//advance(1);
    }
  
  
 
  void InitializeMacroparameters()
  {  const int numMeshes =_meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount(); 
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);  
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[cells]);
	
	//initialize density,velocity  
	for(int c=0; c<nCells;c++)
	  {
	    density[c] =10.0;
	    v[c][0]=10.0;
	    v[c][1]=10.0;
	    v[c][2]=0.0;
	    temperature[c]=5.0;
	    pressure[c]=temperature[c]*density[c];
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
	const int nCells = cells.getCount();
	
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);
	VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[cells]);
	TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[cells]);
	const int N123 = _quadrature.getDirCount(); 
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const TArray& dcxyz= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	
	
	//initialize density,velocity,temperature to zero    
	for(int c=0; c<nCells;c++)
	  {
	    density[c]=0.0;
	    v[c][0]=0.0;
	    v[c][1]=0.0;
	    v[c][2]=0.0;
	    temperature[c]=0.0;   
	  }	
	
	for(int j=0;j<N123;j++){
	  
	  Field& fnd = *_dsfPtr.dsf[j];
	  const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
	  //fprintf(pFile,"%12.6f %E %E \n",dcxyz[j],f[0],density[0]+dcxyz[j]*f[0]);
	  
	  for(int c=0; c<nCells;c++){
	    density[c] = density[c]+dcxyz[j]*f[c];
	    v[c][0]= v[c][0]+(cx[j]*f[c])*dcxyz[j];
	    v[c][1]= v[c][1]+(cy[j]*f[c])*dcxyz[j];
	    v[c][2]= v[c][2]+(cz[j]*f[c])*dcxyz[j];
	    temperature[c]= temperature[c]+(pow(cx[j],2.0)+pow(cy[j],2.0)
					    +pow(cz[j],2.0))*f[c]*dcxyz[j];
	    
	  }
	  
	}
	

	cout << "T" << temperature[0] << endl;
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
	cout << "T" << temperature[0] << endl;
	
      }
    //fclose(pFile);
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
	double Tmuref=273.15;
	double mu_w=0.81; double muref=2.117e-5; //Argon
	double rho_init=1.6034e-4; double T_init= 300; 
	double R=8314.0/40.0;
	double nondim_length=1.0;
	double mu0=rho_init*R*T_init*nondim_length/pow(2*R*T_init,0.5);  
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	
	TArray& density = dynamic_cast<TArray&>(_macroFields.density[cells]);
	TArray& viscosity = dynamic_cast<TArray&>(_macroFields.viscosity[cells]);
	TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[cells]);

	TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[cells]);
	for(int c=0; c<nCells;c++)
	  {
	    viscosity[c]=muref*pow(temperature[c]/Tmuref,mu_w); // viscosity power law
	    collisionFrequency[c]=density[c]*temperature[c]/(muref/mu0*pow(temperature[c]/Tmuref*T_init,mu_w))  ;
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
	double pi(3.14159);
	
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	cout << "temperature [0] = " << temperature[0] << endl;
	cout << "v [0][0] = " << v[0][0] << endl;
	cout << "v [0][1] = " << v[0][1] << endl;
	cout << "v [0][2] = " << v[0][2] << endl;
	for(int j=0;j< numFields;j++){
	  Field& fnd = *_dsfEqPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=density[c]/pow((pi*temperature[c]),1.5)*
	      exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
		    pow((cz[j]-v[c][2]),2.0))/temperature[c]);
	  } 
	  
	}
      }
  }

  void weightedMaxwellian(double weight1,double vel1,double vel2)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getCount();
	double pi(acos(-1.0));
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	for(int j=0;j< numFields;j++){
	  Field& fnd = *_dsfPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=weight1*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-vel1),2.0)+pow((cy[j]-0.0),2.0)+pow((cz[j]-0.0),2.0))/1.0)
	      +(1-weight1)*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-vel2),2.0)+pow((cy[j]-0.0),2.0)+pow((cz[j]-0.0),2.0))/1.0);
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
  
  FlowBCMap& getBCMap() 
    {
      return _bcMap;
    }
  FlowVCMap& getVCMap()
    {
return _vcMap;
    }
  
  KineticModelOptions<T>& 
    getOptions() 
    {
      return _options;
    }
  
  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
	
        const FlowVC<T>& vc = *_vcMap[mesh.getID()];
        
        const StorageSite& cells = mesh.getCells();
	
        shared_ptr<VectorT3Array> vCell(new VectorT3Array(cells.getCount()));

        VectorT3 initialVelocity;
        initialVelocity[0] = _options["initialXVelocity"];
        initialVelocity[1] = _options["initialYVelocity"];
        initialVelocity[2] = _options["initialZVelocity"];
        *vCell = initialVelocity;
        _macroFields.velocity.addArray(cells,vCell);

        
        shared_ptr<TArray> pCell(new TArray(cells.getCount()));
        *pCell = _options["operatingPressure"];
        _macroFields.pressure.addArray(cells,pCell);
     

        shared_ptr<TArray> rhoCell(new TArray(cells.getCount()));
        *rhoCell = vc["density"];
        _macroFields.density.addArray(cells,rhoCell);

        shared_ptr<TArray> muCell(new TArray(cells.getCount()));
        *muCell = vc["viscosity"];
        _macroFields.viscosity.addArray(cells,muCell);

        shared_ptr<TArray> tempCell(new TArray(cells.getCount()));
        *tempCell = _options["operatingTemperature"];
	_macroFields.temperature.addArray(cells,tempCell);
	
	shared_ptr<TArray> collFreqCell(new TArray(cells.getCount()));
        *collFreqCell = vc["viscosity"];
	_macroFields.collisionFrequency.addArray(cells,collFreqCell);

      }
    _niters  =0;
    _initialKmodelNorm = MFRPtr();
  
  }
  
  
  
  void SetBoundaryConditions()
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
		      bc->bcType = "WallBC";
		  }
		else if ((fg.groupType == "velocity-inlet"))
		  {
		    bc->bcType = "WallBC";
		  }
		else if ((fg.groupType == "pressure-inlet") ||
			 (fg.groupType == "pressure-outlet"))
		  {
		    bc->bcType = "PressureBC";
		  }
		else if ((fg.groupType == "symmetry"))
		  {
		    bc->bcType = "SymmetryBC";
		    }
		else if((fg.groupType =="copy "))
		  {
		      bc->bcType = "CopyBC";
		  }
		else
		  throw CException("FlowModel: unknown face group type "
				     + fg.groupType);
	      }
	  }
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
  } 
  
  
  void linearizeKineticModel(LinearSystem& ls, int direction)
  {
    // _velocityGradientModel.compute();
    DiscrList discretizations;
    
    //shared_ptr<Discretization> cd(new ConvectionDiscretization_Kmodel<VectorT3,DiagTensorT3,T> (_meshes,_geomFields,
    //     _flowFields.velocity, _flowFields.massFlux,_flowFields.continuityResidual, _flowFields.velocityGradient));
    //discretizations.push_back(cd);
    //shared_ptr<Discretization> ud(new Underrelaxer<VectorT3,DiagTensorT3,T>
    //     (_meshes,_flowFields.velocity, _options["momentumURF"]));
    //discretizations.push_back(ud);
    
    Field& fnd = *_dsfPtr.dsf[direction]; 
    Field& feq = *_dsfEqPtr.dsf[direction]; 

    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);

    //const double cx1 = dynamic_cast<const double>(cx[direction]);
    //const double cy1 = dynamic_cast<const double>(cy[direction]);
    //const double cz1 = dynamic_cast<const double>(cz[direction]);

    shared_ptr<Discretization>
      sd(new CollisionTermDiscretization<T,T,T>
	 (_meshes, _geomFields, 
	  fnd,feq,
	  _macroFields.collisionFrequency));
    discretizations.push_back(sd);

    /*    shared_ptr<Discretization> 
      cd(new ConvectionDiscretization_Kmodel<T,T,T> 
	 (_meshes,
	  _geomFields,
	  fnd,
	  cx[direction],
	  cy[direction],
 	  cz[direction]));
    discretizations.push_back(cd);
    */
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
	      _options["nonDimLength"],
	      _options.timeDiscretizationOrder));
	
	discretizations.push_back(td);
	
      }
   
    Linearizer linearizer;
    
    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
			 ls.getX(), ls.getB());
  }
  
  void updateTime()
  {
    const int numMeshes = _meshes.size();
    for (int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh = *_meshes[n];

	const StorageSite& cells = mesh.getCells();
	const int numFields= _quadrature.getDirCount(); 

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
      }
  }	  
  
  
  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)  
      {
	const int N123 =_quadrature.getDirCount();

	MFRPtr rNorm;
	for(int direction=0; direction<N123;direction++)
	  {
	    LinearSystem ls;
	    initKineticModelLinearization(ls, direction);
	    ls.initAssembly();
	    linearizeKineticModel(ls,direction);
	    ls.initSolve();
	    
	    
	    MFRPtr iNorm(_options.getKineticLinearSolver().solve(ls));
	    //if (rNorm)
	    //  *rNorm += *iNorm;
	    //else
	    //  rNorm = iNorm;
	     ls.postSolve();
	     ls.updateSolution();
	     _options.getKineticLinearSolver().cleanup();
	  }
	
	//ComputeMacroparameters();
	//ComputeCollisionfrequency();	
	//	_dsfEqPtr.initializeMaxwellian(_macroFields,_dsfEqPtr);
	//	initializeMaxwellianEq();
	  	/*
	if (!_initialKmodelNorm) _initialKmodelNorm = rNorm; 
	MFRPtr normRatio((*rNorm)/(*_initialKModelNorm));
	      cout << _niters << ": " << *rNorm << endl;   
	      _options.getKineticLinearSolver().cleanup();
	      
	      //
	      _niters++;
	      if (*rNorm < _options.absoluteTolerance ||
	      *normRatio < _options.relativeTolerance)
	      break;
	*/	    
	  
      }
    //char * filename="f.txt";
    //itoa(n,filename,10);
    //_dsfPtr.OutputDsf(_dsfPtr,filename);
  }
  
void OutputDsfBLOCK(const char* filename)
  {
    FILE * pFile;
    pFile = fopen(filename,"w"); 
    int N1=_quadrature.getNVCount();
    int N2=_quadrature.getNthetaCount();
    int N3=_quadrature.getNphiCount();
    fprintf(pFile,"%s \n", "VARIABLES= cx, cy, cz, f,fEq,");
    fprintf(pFile, "%s %i %s %i %s %i %s \n","ZONE I=", N3,",J=",N2,",K=",N1,
    	    ",F=BLOCK, VARLOCATION=(NODAL,NODAL,NODAL,NODAL,NODAL)");
    const int numMeshes = _meshes.size();
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
	fprintf(pFile,"%E\n",f[0]);
      } 
      for(int j=0;j< numFields;j++){
	Field& fEqnd = *_dsfEqPtr.dsf[j];
	TArray& fEq = dynamic_cast< TArray&>(fEqnd[cells]);
	fprintf(pFile,"%E\n",fEq[0]);
      }
    }
    fclose(pFile);
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
    for (int n=0; n<numMeshes; n++){
      const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
      const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
      const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
      const int numFields= _quadrature.getDirCount(); 	
      
      const Mesh& mesh = *_meshes[n];
      const StorageSite& cells = mesh.getCells();
      for(int j=0;j< numFields;j++){
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& f = dynamic_cast< TArray&>(fnd[cells]);	
	Field& fEqnd = *_dsfEqPtr.dsf[j];
	TArray& fEq = dynamic_cast< TArray&>(fEqnd[cells]);
	fprintf(pFile,"%E %E %E %E %E\n",cx[j],cy[j],cz[j],f[0],fEq[0]);
      }
    }
  }

  
 private:
  //shared_ptr<Impl> _impl;
 
  const GeomFields& _geomFields;
  const Quadrature<T>& _quadrature;
 
  MacroFields& _macroFields;
  DistFunctFields<T> _dsfPtr;  
  DistFunctFields<T> _dsfPtr1;
  DistFunctFields<T> _dsfPtr2;
  DistFunctFields<T> _dsfEqPtr;
  FlowBCMap _bcMap;
  FlowVCMap _vcMap;

  KineticModelOptions<T> _options;
  int _niters;

  MFRPtr _initialKmodelNorm;
  shared_ptr<Field> _previousVelocity;
  shared_ptr<Field> _KmodelApField;
};

#endif
