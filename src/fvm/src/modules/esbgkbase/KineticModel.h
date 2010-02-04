#ifndef _KINETICMODEL_H_
#define _KINETICMODEL_H_

#include <math.h>
#include "Model.h"
#include "Array.h"
#include "Vector.h"
#include "quadrature.h"

template<class T>
class KineticModel : public Model
{
 public:
  typedef Vector<T,3> VectorT3; 
  
  /**
    * Calculation of macro-parameters density, temperature, components of velocity, pressure
    * by taking moments of distribution function using quadrature points and weights from quadrature.h
    */

  void ComputeMacroparameters(const Mesh& mesh,const MacroParameters& macroPr, 
			      const Quadrature<T>& quad, const DistFunctFields& distfunfields) 
  {
    const StorageSite& cells = mesh.getCells();
    const int nCells = cells.getSelfCount();
    double pi(3.14159);
    
    const TArray& density = dynamic_cast<const TArray&>(macroPr.density[cells]);
    const TArray& temperature = dynamic_cast<const TArray&>(macroPr.temperature[cells]);
    const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(macroPr.velocity[cells]);
    
    const int N123 = quad.getDirCount(); 
    const TArray& cx = dynamic_cast<const TArray&>(quad.cxPtr[cells]);
    const TArray& cy = dynamic_cast<const TArray&>(quad.cyPtr[cells]);
    const TArray& cz = dynamic_cast<const TArray&>(quad.czPtr[cells]);
    const TArray& dcxyz= dynamic_cast<const TArray&>(quad.dcxyzPtr[cells]);
    
    double loc1,loc2,loc3,loc4;
    for(int c=0; c<ncells;c++){
     
      loc1=0.0;loc2=0.0;loc3=0.0;
      for(int j=0;j<N123;j++){
	loc1=loc1+f[j]*dcxyz[j];
	loc2=loc2+(cx[j]*f[j])*dcxyz[j];
	loc3=loc3+(cy[j]*f[j])*dcxyz[j];
	loc4=loc4+(cz[j]*f[j])*dcxyz[j];
      }
      
      density[c]=loc1;
      v[c][0]=loc2/density[c];
      v[c][1]=loc3/density[c];
      v[c][2]=loc4/density[c];
      
      loc2=0;loc3=0; 
      for(j=0;j<N123;j++){
	loc2=loc2+(pow(cx[j]-v[c][0],2.0)+pow(cy[j]-v[c][1],2.0)+pow(cz[j]-v[c][2],2.0))*f[j]*dcxyz[j];
	loc3=loc3+pow(cx[j]-v[c][0]],2.0)*f[j]*dcxyz[j];
    }

    //loop through cells 
    temperature[c]=loc2/(1.5*density[c]);                            
    //xTemp[c]] =loc3/(1.5*density[c]);
  }
  
  
  /** 
   * Collision frequency
   */
  void ComputeCollisionfrequency(Const Mesh& mesh)
  {
    
  }
  
 private:
  shared_ptr<Impl> _impl;
};

#endif

