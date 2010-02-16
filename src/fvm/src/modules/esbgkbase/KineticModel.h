#ifndef _KINETICMODEL_H_
#define _KINETICMODEL_H_

#include <math.h>
#include "Model.h"
#include "Array.h"
#include "Vector.h"
#include "quadrature.h"
#include "DistFunctFields.h"
#include "MacroParameters.h"
#include "FlowFields.h"

template<class T>
class KineticModel : public Model
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3; 
  typedef Array<VectorT3> VectorT3Array;
  
  /**
   * Calculation of macro-parameters density, temperature, components of velocity, pressure
   * by taking moments of distribution function using quadrature points and weights from quadrature.h
   */
  DistFunctFields<T>& distfunfields;
  MacroParameters& macroPr;
  KineticModel(const Mesh& mesh,const Quadrature<T>& quad)
    {
      //_mesh(mesh),
      //_quad(quad),
      //const DistFunctFields<T>& distfunfields,
     DistFunctFields<T>&  distfuncfields(const Mesh &mesh, const MacroParameters& macroPr,const Quadrature<T>& quad);   
   
    }
  
  void ComputeMacroparameters(const Mesh& mesh,const MacroParameters& macroPr, 
			      const Quadrature<T>& quad, const DistFunctFields<T>& distfunfields) 
  {
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();
    double pi(3.14159);
    
    const TArray& density = dynamic_cast<const TArray&>(macroPr.density[cells]);
    const TArray& temperature = dynamic_cast<const TArray&>(macroPr.temperature[cells]);
    const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(macroPr.velocity[cells]);
     const TArray& pressure = dynamic_cast<const TArray&>(macroPr.pressure[cells]);
    const int N123 = quad.getDirCount(); 
    const TArray& cx = dynamic_cast<const TArray&>(quad.cxPtr[cells]);
    const TArray& cy = dynamic_cast<const TArray&>(quad.cyPtr[cells]);
    const TArray& cz = dynamic_cast<const TArray&>(quad.czPtr[cells]);
    const TArray& dcxyz= dynamic_cast<const TArray&>(quad.dcxyzPtr[cells]);
    
    //initialize density,velocity,temperature to zero    
    for(int c=0; c<nCells;c++)
      {
	density[c] =0.0;
	v[c][0]=0.0;
	v[c][1]=0.0;
	v[c][2]=0.0;
	temperature[c]=0.0;
	
      }	
    
    //loop through directions
    for(int j=0;j<N123;j++){
      const TArray& f = dynamic_cast<const TArray&>(distfunfields.dsf[j]);
      //loop through cells
      for(int c=0; c<nCells;c++)
	{
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
  
  /**
   * Collision frequency
   */
  void ComputeCollisionfrequency(const Mesh& mesh,const MacroParameters& macroPr)  {
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();
    double pi(3.14159);
    
    const TArray& viscosity = dynamic_cast<const TArray&>(macroPr.viscosity[cells]);
    //viscosity[c]=muref*pow(temperature[c]/Tmuref,mu_w); // viscosity power law
    Field CollisionFrequency;
    
  }

  void MacroParametersToFlowfields(const Mesh& mesh,const MacroParameters& macroPr,
 const FlowFields& flowfields, int method)
  {
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();
    const TArray& d_macro = dynamic_cast<const TArray&>(macroPr.density[cells]);
    const TArray& d_ffields = dynamic_cast<const TArray&>(flowfields.density[cells]);
    const TArray& t_macro = dynamic_cast<const TArray&>(macroPr.density[cells]);
    const TArray& t_ffields = dynamic_cast<const TArray&>(flowfields.density[cells]);
    const VectorT3Array& v_macro = dynamic_cast<const VectorT3Array&>(macroPr.velocity[cells]);
    const VectorT3Array& v_ffields = dynamic_cast<const VectorT3Array&>(flowfields.velocity[cells]);
    if(method=1){ //copy from Flowfields to  MacroParameters
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
  
