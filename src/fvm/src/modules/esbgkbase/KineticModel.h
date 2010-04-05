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
  
 KineticModel(const MeshList& meshes, FlowFields& ffields,MacroFields& macroFields, const Quadrature<T>& quad):
  //KineticModel(const MeshList& meshes, FlowFields& ffields, const Quadrature<T>& quad):
  
Model(meshes),
  _meshes(meshes),
    _quadrature(quad),
    _macroFields(macroFields),
  _dsfPtr(_meshes,_quadrature)
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
      //ComputeMacroparameters(_meshes,ffields,_quadrature,_dsfPtr);
      ComputeMacroparameters(_meshes,_macroFields,_quadrature,_dsfPtr); //calculate density,velocity,temperature
      ComputeCollisionfrequency(_meshes, _macroFields); //calculate viscosity, collisionFrequency
    }
  

  
 void InitializeMacroparameters(const MeshList& meshes, MacroFields& macroPr,
			     const Quadrature<T>& quad,const DistFunctFields<T>& dsff)
 {  const int numMeshes =meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
        const Mesh& mesh = *meshes[n];
   const StorageSite& cells = mesh.getCells();
   const int nCells = cells.getSelfCount(); 
   
   TArray& density = dynamic_cast<TArray&>(macroPr.density[cells]);  
   VectorT3Array& v = dynamic_cast<VectorT3Array&>(macroPr.velocity[cells]);
   TArray& temperature = dynamic_cast<TArray&>(macroPr.temperature[cells]);
   TArray& pressure = dynamic_cast<TArray&>(macroPr.pressure[cells]);
 
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
 
 void ComputeMacroparameters(const MeshList& meshes, MacroFields& macroPr, 
			  const Quadrature<T>& quad, const DistFunctFields<T>& dsff) 
 {  
   //FILE * pFile;
   //pFile = fopen("distfun_mf.txt","w");
   const int numMeshes = meshes.size();
   for (int n=0; n<numMeshes; n++)
     {
     
       const Mesh& mesh = *meshes[n];
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
	   density[c]=0.0;
	   v[c][0]=0.0;
	   v[c][1]=0.0;
	   v[c][2]=0.0;
	   temperature[c]=0.0;   
	 }	
       
       for(int j=0;j<N123;j++){
	 
	 Field& fnd = *dsff.dsf[j];
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
	 
	 
       }
       
       
     }
   //fclose(pFile);
 }
 
 /*
  * Collision frequency
  *
  *
  */
  void ComputeCollisionfrequency(const MeshList& meshes, MacroFields& macroPr)  {
    const int numMeshes = meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	double Tmuref=273.15;
	double mu_w=0.81; double muref=2.117e-5; //Argon
	double rho_init=1.6034e-4; double T_init= 300; 
	double R=8314.0/40.0;
	double nondim_length=1.0;
	double mu0=rho_init*R*T_init*nondim_length/pow(2*R*T_init,0.5);  
	const Mesh& mesh = *meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();
	
	TArray& density = dynamic_cast<TArray&>(macroPr.density[cells]);
	TArray& viscosity = dynamic_cast<TArray&>(macroPr.viscosity[cells]);
	TArray& temperature = dynamic_cast<TArray&>(macroPr.temperature[cells]);

	TArray& collisionFrequency = dynamic_cast<TArray&>(macroPr.collisionFrequency[cells]);
	for(int c=0; c<nCells;c++)
	  {
	    viscosity[c]=muref*pow(temperature[c]/Tmuref,mu_w); // viscosity power law
	    collisionFrequency[c]=density[c]*temperature[c]/(muref/mu0*pow(temperature[c]/Tmuref*T_init,mu_w))  ;
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
FlowModelOptions<T>& 
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
  

 private:
  //shared_ptr<Impl> _impl;
  const MeshList& _meshes;
  const Quadrature<T>& _quadrature;
 
  MacroFields& _macroFields;
  DistFunctFields<T> _dsfPtr; 
  FlowBCMap _bcMap;
  FlowVCMap _vcMap;

  FlowModelOptions<T> _options;

};

#endif
  
