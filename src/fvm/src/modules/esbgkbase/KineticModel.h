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
#include "GenericKineticBCS.h"

#include "Linearizer.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"

#include "MatrixOperation.h"
#include "NumType.h"

template<class T>
class KineticModel : public Model
{
 public:
  typedef typename NumTypeTraits<T>:: T_Scalar T_Scalar;
  typedef Array<T> TArray;
  typedef  Array2D<T> TArray2D;
  typedef Vector<T,3> VectorT3; 
  typedef Array<VectorT3> VectorT3Array;
  typedef  std::vector<Field*> stdVectorField;
  typedef  DistFunctFields<T> TDistFF;

  typedef Vector<T,5> VectorT5; 
  typedef Array<VectorT5> VectorT5Array;
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
    //KineticModel(const MeshList& meshes, FlowFields& ffields, const Quadrature<T>& quad):
    
    Model(meshes),
    _geomFields(geomFields),
    _quadrature(quad),
    _macroFields(macroFields),
    _dsfPtr(_meshes,_quadrature),
    _dsfPtr1(_meshes,_quadrature),
    _dsfPtr2(_meshes,_quadrature),
    _dsfEqPtr(_meshes,_quadrature),
    _dsfEqPtrES(_meshes,_quadrature)
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
      initializeMaxwellianEq(); //equilibrium distribution
      
      //EquilibriumDistributionBGK();
      callBoundaryConditions();
     
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
	TArray& Txz = dynamic_cast<TArray&>(_macroFields.Txz[cells]);
	TArray& Tyz = dynamic_cast<TArray&>(_macroFields.Tyz[cells]);
	//if ( MPI::COMM_WORLD.Get_rank() == 0 ) {cout << "ncells="<<nCells<<endl;}
	
	//initialize density,velocity  
	for(int c=0; c<nCells;c++)
	  {
	    density[c] =1.0;
	    v[c][0]=0.0;
	    v[c][1]=0.0;
	    v[c][2]=0.0;
	    if(_options.fgamma>0){
	    //BGK
	    coeff[c][0]=1.0/pow(pi,1.5);
	    coeff[c][1]=1.0;
	    coeff[c][2]=0.0;coeff[c][3]=0.0;coeff[c][4]=0.0;
	    }
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
	   

	    temperature[c]=1.0;
	    pressure[c]=temperature[c]*density[c];
	    
	    Txx[c]=0.5;
	    Tyy[c]=0.5;
	    Tzz[c]=0.5;
	    Txy[c]=0.0;
	    Txz[c]=0.0;
	    Tyz[c]=0.0;
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
	const int nCells = cells.getCount();  //only interior
	
	
	TArray& Entropy = dynamic_cast<TArray&>(_macroFields.Entropy[cells]);
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
	    Entropy[c]=0.0;
	  }	
	const T rho_init=_options["rho_init"]; 
const T T_init=_options["T_init"]; 
const T R=8314.0/_options["molecularWeight"];
T u_init=pow(2.0*R*T_init,0.5);
	for(int j=0;j<N123;j++){
	  
	  Field& fnd = *_dsfPtr.dsf[j];
	  const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
	  //fprintf(pFile,"%12.6f %E %E \n",dcxyz[j],f[0],density[0]+dcxyz[j]*f[0]);
	  
	  for(int c=0; c<nCells;c++){
	    density[c] = density[c]+dcxyz[j]*f[c];
	    v[c][0]= v[c][0]+(cx[j]*f[c])*dcxyz[j];
	    v[c][1]= v[c][1]+(cy[j]*f[c])*dcxyz[j];
	    v[c][2]= v[c][2]+(cz[j]*f[c])*dcxyz[j];
	   Entropy[c]=Entropy[c]+f[c]*dcxyz[j]*(1-log(pow(6.26068E-34,3.0)*f[c]/pow(40.0E-26/6.023,4)*rho_init/pow(u_init,3.0)));

		 temperature[c]= temperature[c]+(pow(cx[j],2.0)+pow(cy[j],2.0)
					    +pow(cz[j],2.0))*f[c]*dcxyz[j];
	    
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

	
	
      }
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
	TArray& Txz = dynamic_cast<TArray&>(_macroFields.Txz[cells]);
	TArray& Tyz = dynamic_cast<TArray&>(_macroFields.Tyz[cells]);
	
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
	    Txz[c]=0.0;
	    Tyz[c]=0.0;
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
	    Txy[c]=Txy[c]+(cx[j])*(cy[j])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];
	    Txz[c]=Txz[c]+(cx[j])*(cz[j])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];
	    Tyz[c]=Tyz[c]+(cy[j])*(cz[j])*((1-1/Pr)*f[c]+1/Pr*fgam[c])*wts[j];  
	  }
	}
		
	for(int c=0; c<nCells;c++){
	  Txx[c]=Txx[c]/density[c];
	  Tyy[c]=Tyy[c]/density[c];
	  Tzz[c]=Tzz[c]/density[c];
	  Txy[c]=Txy[c]/density[c];
	  Txz[c]=Txz[c]/density[c];
	  Tyz[c]=Tyz[c]/density[c];
	}
	//if (n==0){cout << "ESBGK-" <<"Txx "<< Txx[0]<<" Tyy "<<Tyy[0]<<" Tzz "<<Tzz[0]<<endl;}
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
	const T nondim_length=_options["nonDimLength"];
	/*
	const T rho_init=9.98e-7; 
	const T T_init= 273.15;
	const T mu_w=0.81;
	const T Tmuref=273.15;
	const T muref=2.117e-5;
	const T R=8314.0/40.0;
	const T nondim_length=1.0;
	*/
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
	    collisionFrequency[c]=density[c]*temperature[c]/viscosity[c]*mu0;//(muref/mu0*pow(temperature[c]/ Tmuref*T_init, mu_w))  ;
	  }
	if(_options.fgamma==2){
	  for(int c=0; c<nCells;c++)
	    {collisionFrequency[c]=_options.Prandtl*collisionFrequency[c];}
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
	//double pi(3.14159);
	const double pi=_options.pi;
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 

	for(int j=0;j< numFields;j++){
	  Field& fndEq = *_dsfEqPtr.dsf[j];
	  TArray& fEq = dynamic_cast< TArray&>(fndEq[cells]);
	  for(int c=0; c<nCells;c++){
	    fEq[c]=density[c]/pow((pi*temperature[c]),1.5)*
	      exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
		    pow((cz[j]-v[c][2]),2.0))/temperature[c]);
	
	  } 
	  
	} 
	
	/*
	if(_options.fgamma==2){
	  for(int j=0;j< numFields;j++){
	    Field& fndEqES = *_dsfEqPtrES.dsf[j];
	    TArray& fEqES = dynamic_cast< TArray&>(fndEqES[cells]);
	    for(int c=0; c<nCells;c++){
	      fEqES[c]=density[c]/pow((pi*temperature[c]),1.5)*
		exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
		      pow((cz[j]-v[c][2]),2.0))/temperature[c]);
	    } 
	    
	    }
	 }
	*/
	
      }
  }
  
  void NewtonsMethodBGK(const int ktrial)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {		
	//cout << "inside NewtonsMethod" <<endl;
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
	    // TArray* fvecPtr;
	    //fvecPtr=new TArray(5);
	    //TArray & fvec= *fvecPtr;
	    
	    fvec[0]=density[c];
	    fvec[1]=density[c]*v[c][0];
	    fvec[2]=density[c]*v[c][1];
	    fvec[3]=density[c]*v[c][2];
	    fvec[4]=1.5*density[c]*temperature[c]+density[c]*(pow(v[c][0],2)+pow(v[c][1],2)+pow(v[c][2],2.0));
	    //calculate Jacobian
	    
	    setJacobianBGK(fjac,fvec,coeff[c],v[c],c);
	    
	    
	    /*
	    if(c==0){
	      cout << "ktrial"<< trial <<endl;
	      cout << "fj  " <<fjac(1,1) <<endl;
	      cout << "fv  " <<fvec[1] <<endl;}
	    */
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
	     
	    //if(c==20){cout<<"fjac " <<fjac(1,1) << endl;}
	    
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
	    
	    /*
	    if (c==_options.printCellNumber) { cout <<"trial" <<trial <<endl;
	      cout << "BGK-CELL " <<c << " c0 "<< coeff[c][0]<<"  c1 " <<coeff[c][1]<< " c2 "<< coeff[c][2]<< " c3 "<< coeff[c][3] << " c4 "<< coeff[c][4]<<endl;}
	    */
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
	//	const double pi(3.14159);
	const double pi=_options.pi;
	const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[cells]);
	const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[cells]);
	
	//initialize coeff
	VectorT5Array& coeff = dynamic_cast<VectorT5Array&>(_macroFields.coeff[cells]);
	for(int c=0; c<nCells;c++){
	  coeff[c][0]=density[c]/pow((pi*temperature[c]),1.5);
	  coeff[c][1]=1/temperature[c];
	  coeff[c][2]=0.0;
	  coeff[c][3]=0.0;
	  coeff[c][4]=0.0;	    
	}
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	//call Newtons Method
	//const int ktrial(5);
	//cout << "calling Newtons method" <<endl;
	NewtonsMethodBGK(ktrial);
	//int cellno=_options.printCellNumber;
	//if (n==0){cout << " coeffs-BGK " << coeff[cellno][0] <<" cx^2  " <<coeff[cellno][1] <<" cx  " <<coeff[cellno][2] <<endl;}

	//calculate perturbed maxwellian for BGK
	//const VectorT5Array& coeff = dynamic_cast<const VectorT5Array&>(_macroFields.coeff[cells]);
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
      //if (c==20){cout<< "E,malpha " << Econst<<malphaBGK(j,0)<<"fvec" << fvec[0]<<endl;}
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
	const TArray& Txz = dynamic_cast<const TArray&>(_macroFields.Txz[cells]);
	const TArray& Tyz = dynamic_cast<const TArray&>(_macroFields.Tyz[cells]);
	
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
    	
	VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
	
	for(int c=0; c<nCells;c++){
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
	    fvec[8]=density[c]*(v[c][0]*v[c][2]+Txz[c]);
	    fvec[9]=density[c]*(v[c][1]*v[c][2]+Tyz[c]);
	    
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


	const VectorT5Array& coeff = dynamic_cast<const VectorT5Array&>(_macroFields.coeff[cells]);
	//initialize coeffg
	VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[cells]);
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
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[cells]);
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	
	//cout <<"calling Newtons Method"<<endl;//call Newtons Method
	//const int ktrial(1);
	NewtonsMethodESBGK(ktrial);
	//cout <<"called Newtons Method"<<endl;

	//calculate perturbed maxwellian for BGK
	
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
				    +coeffg[c][7]*cx[j]*cy[j]+coeffg[c][8]*cx[j]*cz[j]
				    +coeffg[c][9]*cy[j]*cz[j]);
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
			 +xn[7]*cx[j]*cy[j]+xn[8]*cx[j]*cz[j]+xn[9]*cy[j]*cz[j])*wts[j];     
      
      for (int row=0;row<10;row++){
	fvec[row]+= -Econst*malphaESBGK(j,row); //smm
	//	if (c==0 && j<10){cout <<"dir row"<< j << row << fvec[row] <<endl;}
      }
      
      mexp[0]=-Econst/xn[0];
      mexp[1]=Econst*pow(Cc1,2);
      mexp[2]=-Econst*Cc1;
      mexp[3]=Econst*pow(Cc2,2);
      mexp[4]=-Econst*Cc2; 
      mexp[5]=Econst*pow(Cc3,2);
      mexp[6]=-Econst*Cc3;
      
      mexp[7]=-Econst*cx[j]*cy[j];
      mexp[8]=-Econst*cx[j]*cz[j];
      mexp[9]=-Econst*cy[j]*cz[j];
   
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
	//double pi(3.14159);
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

  void weightedMaxwellian(double weight1,double vel1,double vel2,double temp1,double temp2)
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
	    f[c]=weight1*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-vel1),2.0)+pow((cy[j]-0.0),2.0)+pow((cz[j]-0.0),2.0))/temp1)
	      +(1-weight1)*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-vel2),2.0)+pow((cy[j]-0.0),2.0)+pow((cz[j]-0.0),2.0))/temp2);
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
  
  void init()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *_meshes[n];
	
        const KineticVC<T>& vc = *_vcMap[mesh.getID()];
        
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
	

	//coeffs for perturbed BGK distribution function
	shared_ptr<VectorT5Array> coeffCell(new VectorT5Array(cells.getCount()));
        VectorT5 initialCoeff;
        initialCoeff[0] = 1.0;//vc["density"]/pow(_options["pi"]*_options["operatingTemperature"],1.5);
        initialCoeff[1] = 1.0;//1/_options["operatingTemperature"];
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

	shared_ptr<TArray> tempxzCell(new TArray(cells.getCount()));
        *tempxzCell = 0.0;
	_macroFields.Txz.addArray(cells,tempxzCell);

	shared_ptr<TArray> tempyzCell(new TArray(cells.getCount()));
        *tempyzCell = 0.0;
	_macroFields.Tyz.addArray(cells,tempyzCell);

	

        shared_ptr<TArray> EntropyCell(new TArray(cells.getCount()));
        *EntropyCell = 0.0;
        _macroFields.Entropy.addArray(cells,EntropyCell);
      }
    _niters  =0;
    _initialKmodelNorm = MFRPtr();
    //_initialKmodelvNorm = MFRPtr();
  
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
		else if (fg.groupType == "velocity-inlet")
		  {
		    bc->bcType = "InletBC";
		  }
		else if (fg.groupType == "pressure-inlet")
		  {
		    bc->bcType = "PressureInletBC";
		  }
		else if (fg.groupType == "pressure-outlet")
		  {
		    bc->bcType = "PressureOutletBC";
		  }
		else if (fg.groupType == "velocity-inlet")
		  {
		    bc->bcType = "VelocityInletBC";
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
 	  cz[direction]));
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
	      _options["nonDimLength"],
	      _options.timeDiscretizationOrder));
	
	discretizations.push_back(td);
	
      }
   
    Linearizer linearizer;
    
    linearizer.linearize(discretizations,_meshes,ls.getMatrix(),
			 ls.getX(), ls.getB());

    // boundary conditions
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
		  {
		    gkbc.applyExtrapolationBC(f);
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
		  if(c_dot_en -wallV_dot_en < T_Scalar(0.0))
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
	callBoundaryConditions();  //new
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

	ComputeMacroparameters();	//update macroparameters
        ComputeCollisionfrequency();
	if (_options.fgamma==0){initializeMaxwellianEq();}
	else{ EquilibriumDistributionBGK();}
	
	if (_options.fgamma==2){EquilibriumDistributionESBGK();}
	
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
	    
	    FloatValEvaluator<VectorT3>
	      bVelocity(bc.getVal("specifiedXVelocity"),
			bc.getVal("specifiedYVelocity"),
			bc.getVal("specifiedZVelocity"),
			faces);
	    FloatValEvaluator<T>
	      bTemperature(bc.getVal("specifiedTemperature"),
			   faces);
	    FloatValEvaluator<T>
	      bPressure(bc.getVal("specifiedPressure"),
			   faces);

	    if (bc.bcType == "WallBC")
	      {	
		//cout <<"applying wallbc"<<endl;
		kbc.applyDiffuseWallBC(bVelocity,bTemperature);
		//cout <<"applied"<<endl;
	      }
	    else if(bc.bcType=="SymmetryBC")
	      {
		//cout <<"applying specularbc"<<endl;
		kbc.applySpecularWallBC();
		//cout <<"applied"<<endl;
	      } 
	    else if(bc.bcType=="ZeroGradBC")
	      {
		kbc.applyZeroGradientBC();
	      }
	    else if(bc.bcType=="PressureInletBC")
	      {
		//cout <<"applying pressureinlet"<<endl;
		kbc.applyPressureInletBC(bTemperature,bPressure);
		//cout <<"applied"<<endl;
	     } 
	    else if(bc.bcType=="VelocityInletBC")
	       {
		kbc.applyVelocityInletBC(bTemperature,bVelocity);
	     }
	    else if(bc.bcType=="PressureOutletBC")
	      {
	    	kbc.applyPressureOutletBC(bTemperature,bPressure);
	      }
	    
	  }
      }
   }

 
  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)  
      {
	const int N123 =_quadrature.getDirCount();
	MFRPtr rNorm;
	//MFRPtr vNorm;
       	//const TArray& cx= dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	//const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

	//callBoundaryConditions();
   	  _macroFields.velocity.syncLocal();
 	  _macroFields.temperature.syncLocal();
  	  _macroFields.density .syncLocal();

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

	
	
	if (!_initialKmodelNorm) _initialKmodelNorm = rNorm;
	//if (!_initialKmodelvNorm) _initialKmodelvNorm = vNorm;
	if (_niters < 5)
        {
             _initialKmodelNorm->setMax(*rNorm);
	     // _initialKmodelvNorm->setMax(*vNorm);
            
        } 
 
	MFRPtr normRatio((*rNorm)/(*_initialKmodelNorm));	
	//	MFRPtr vnormRatio((*vNorm)/(*_initialKmodelvNorm));
	//if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{cout << _niters << ": " << *rNorm <<endl; }

	_niters++;
	if ((*rNorm < _options.absoluteTolerance)||(*normRatio < _options.relativeTolerance )) //&& ((*vNorm < _options.absoluteTolerance)||(*vnormRatio < _options.relativeTolerance )))
	  break;
	
	callBoundaryConditions();
	ComputeMacroparameters();	//update macroparameters
        ComputeCollisionfrequency();
	//update equilibrium distribution function 0-maxwellian, 1-BGK,2-ESBGK
	if (_options.fgamma==0){initializeMaxwellianEq();}
	else{ EquilibriumDistributionBGK();}
	if (_options.fgamma==2){EquilibriumDistributionESBGK();}	    

      }


    //char * filename="f.txt";
    //itoa(n,filename,10);
    //_dsfPtr.OutputDsf(_dsfPtr,filename);
  }

  void advance_dir(const int niter)
  {
    const int N123 =_quadrature.getDirCount();
   
    MFRPtr rNorm;
	  
    for(int direction=0; direction<N123;direction++)
      {


	for(int n=0; n<niter; n++)  
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
		rArray0 = rArrayd;//*wts[direction];
		
	      }
	    ls.postSolve();
	    ls.updateSolution();
	    
	    _options.getKineticLinearSolver().cleanup();
	      
	    
	    if (!_initialKmodelNorm) _initialKmodelNorm = rNorm;
	    
	    if (_niters < 5)
	      {
		_initialKmodelNorm->setMax(*rNorm);	
	      } 
	    
	    MFRPtr normRatio((*rNorm)/(*_initialKmodelNorm));	
	    
	    if(direction ==_options["printDirectionNumber"]){cout << _niters << ": " << *kNorm <<endl; }
	    
	   
	    _niters++;
	    if ((*rNorm < _options.absoluteTolerance)||(*normRatio < _options.relativeTolerance )) 
	      break;
	  
	    
	  }

	  
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
	fprintf(pFile,"%E\n",f[0]);
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
  int _niters;

  MFRPtr _initialKmodelNorm;
  MFRPtr _initialKmodelvNorm;
  shared_ptr<Field> _previousVelocity;
  shared_ptr<Field> _KmodelApField;
};

#endif
