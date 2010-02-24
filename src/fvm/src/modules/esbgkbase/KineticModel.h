#ifndef _KINETICMODEL_H_
#define _KINETICMODEL_H_

#include <math.h>
#include "Model.h"
#include "Array.h"
#include "Vector.h"
#include "Mesh.h"
#include "quadrature.h"
#include "DistFunctFields.h"
#include "MacroParameters.h"
#include "FlowFields.h"

template<class T>
class KineticModel 
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3; 
  typedef Array<VectorT3> VectorT3Array;
  typedef  std::vector<Field*> stdVectorField;
  typedef  DistFunctFields<T> TDistFF;
  /**
   * Calculation of macro-parameters density, temperature, components of velocity, pressure
   * by taking moments of distribution function using quadrature points and weights from quadrature.h
   */
  TDistFF* dsfPtr;
 //MacroParameters& macroPr;
  KineticModel(const Mesh& mesh, FlowFields& ffields,const Quadrature<T>& quad)
    {
      
      //dsfPtr=new TDistFF(mesh,macroPr,quad);   
    dsfPtr = new TDistFF(mesh, quad);      
    ComputeMacroparameters(mesh,ffields,quad,*dsfPtr);
    }
  
 void ComputeMacroparameters(const Mesh& mesh, FlowFields& macroPr,
			     const Quadrature<T>& quad,const DistFunctFields<T>& dsff)
 {
   const StorageSite& cells = mesh.getCells();
   const int nCells = cells.getSelfCount(); 
   
   TArray& density = dynamic_cast<TArray&>(macroPr.density[cells]);  
   VectorT3Array& v = dynamic_cast<VectorT3Array&>(macroPr.velocity[cells]);
   const int N123 = quad.getDirCount(); 
   
   const TArray& cx = dynamic_cast<const TArray&>(*quad.cxPtr);
   const TArray& cy = dynamic_cast<const TArray&>(*quad.cyPtr);
   const TArray& cz = dynamic_cast<const TArray&>(*quad.czPtr);
   const TArray& dcxyz= dynamic_cast<const TArray&>(*quad.dcxyzPtr);
   
   //initialize density,velocity  
   for(int c=0; c<nCells;c++)
     {
       density[c] =1.0/N123;
       v[c][0]=1.0;
       v[c][1]=1.0;
       v[c][2]=0.0;	
     }	

   for(int j=0;j<N123;j++){  //loop through directions
     Field& fnd = *dsff.dsf[j];
     const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
     
     for(int c=0; c<nCells;c++){  //loop through cells
       	 density[c] =density[c]+f[c]*dcxyz[j];
	 v[c][0]=v[c][0]+(cx[j]*f[c])*dcxyz[j];
	 v[c][1]=v[c][1]+(cy[j]*f[c])*dcxyz[j];
	 v[c][2]=v[c][2]+(cz[j]*f[c])*dcxyz[j];
       }	
   }
   
   for(int c=0; c<nCells;c++){
     v[c][0]=v[c][0]/density[c];
     v[c][1]=v[c][1]/density[c];
     v[c][2]=v[c][2]/density[c];   
     }
   
 }
 
 T ComputeMacroparameters(const Mesh& mesh, MacroParameters& macroPr, 
			  const Quadrature<T>& quad, const DistFunctFields<T>& dsff) 
 {
   const StorageSite& cells = mesh.getCells();
   const int nCells = cells.getSelfCount();
   
   
   TArray& density = dynamic_cast<TArray&>(macroPr.density[cells]);
   TArray& temperature = dynamic_cast<TArray&>(macroPr.temperature[cells]);
   VectorT3Array& v = dynamic_cast<VectorT3Array&>(macroPr.velocity[cells]);
   TArray& pressure = dynamic_cast<TArray&>(macroPr.pressure[cells]);
   const int N123 = quad.getDirCount(); 
   const TArray& cx = dynamic_cast<const TArray&>(*quad.cxPtr);
   const TArray& cy = dynamic_cast<const TArray&>(*quad.cyPtr);
   const TArray& cz = dynamic_cast<const TArray&>(*quad.czPtr);
   const TArray& dcxyz= dynamic_cast<const TArray&>(*quad.dcxyzPtr);
   
   
   //initialize density,velocity,temperature to zero    
   for(int c=0; c<nCells;c++)
     {
       density[c] =0.0;
       v[c][0]=0.0;
       v[c][1]=0.0;
       v[c][2]=0.0;
       temperature[c]=0.0;
       
     }	
 
   for(int j=0;j<N123;j++){
    
     Field& fnd = *dsff.dsf[j];
     const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
     
     for(int c=0; c<nCells;c++){
	 density[c] =density[c]+f[c]*dcxyz[j];
	 v[c][0]=v[c][0]+(cx[j]*f[c])*dcxyz[j];
	 v[c][1]=v[c][1]+(cy[j]*f[c])*dcxyz[j];
	 v[c][2]=v[c][2]+(cz[j]*f[c])*dcxyz[j];
	 temperature[c]=temperature[c]+(pow(cx[j],2.0)+pow(cy[j],2.0)
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
     //viscosity[c]=muref*pow(temperature[c]/Tmuref,mu_w); // viscosity power law
     
     //loc2=loc2+(pow(cx[j]-v[c][0],2.0)+pow(cy[j]-v[c][1],2.0)+pow(cz[j]-v[c][2],2.0))*f[j]*dcxyz[j];
     //loc3=loc3+pow(cx[j]-v[c][0]],2.0)*f[j]*dcxyz[j];
   }
   
   
 }
 
 /*
  * Collision frequency
  *
  *
  void ComputeCollisionfrequency(const Mesh& mesh,MacroParameters& macroPr)  {
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();
    double pi(3.14159);
    
    TArray& viscosity = dynamic_cast<TArray&>(macroPr.viscosity[cells]);
    //viscosity[c]=muref*pow(temperature[c]/Tmuref,mu_w); // viscosity power law
    //Field CollisionFrequency;
    
  }
  */

  
  T MacroParametersToFlowfields(const Mesh& mesh, MacroParameters& macroPr,
  FlowFields& flowfields, int method)
  {
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();
    TArray& d_macro = dynamic_cast<TArray&>(macroPr.density[cells]);
    TArray& d_ffields = dynamic_cast<TArray&>(flowfields.density[cells]);
    TArray& t_macro = dynamic_cast<TArray&>(macroPr.density[cells]);
    TArray& t_ffields = dynamic_cast<TArray&>(flowfields.density[cells]);
    VectorT3Array& v_macro = dynamic_cast<VectorT3Array&>(macroPr.velocity[cells]);
    VectorT3Array& v_ffields = dynamic_cast<VectorT3Array&>(flowfields.velocity[cells]);
    if(method==1){ //copy from Flowfields to  MacroParameters
      for(int c=0; c<nCells;c++){
	d_macro[c]=d_ffields[c]; 
	t_macro[c] =t_ffields[c];
	v_macro[c]=v_ffields[c];}
    }
    else {   //copy from MacroParameters to Flowfields
      for(int c=0; c<nCells;c++){
	d_ffields =	d_macro[c];
	t_ffields[c]=t_macro[c];
	v_ffields[c]=v_macro[c];}
    }
  }
  
  
 private:
  //shared_ptr<Impl> _impl;
};

#endif
  
