#ifndef _DISTFUNCTFIELDS_H_
#define _DISTFUNCTFIELDS_H_

#include "Mesh.h"
#include "Quadrature.h"
#include <stdio.h>
#include "Vector.h"

#include "Field.h"
#include "MacroFields.h"
#include "FlowFields.h"

#include "FlowBC.h"

#include <math.h>

template <class T>

/**
 * Class DistFunctFields for Distribution function in ESBGK simulations
 * A collection of distribution function fields in all N123 directions 
 * velocity quadrature points defined in quadrature.h
 */
class DistFunctFields
{
 public:
  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  
  //std::vector<shared_ptr<Field> > DistFunctFieldsPtr;
  std::vector<Field*> dsf;
 

 /**
   * Maxwellian Distribution - sets one direction as one field
   * ix, iy - spatial location
   * j - discrete velocity ordinate
   * f[ix][iy][j]=Rho[ix][iy]/pow((pi*Temp[ix][iy]),1.5)*
   *             exp(-((pow(cx[j]-xVel[ix][iy]),2.0)
   *            +pow((cy[j]-yVel[ix][iy]),2.0)+pow(cz[j],2.0))/Temp[ix][iy])
   * @param mesh - spatial mesh with cells,faces
   * @param macroPr - Fields of macroparameters such as density, temperature, velocity, pressure, viscosity
   * @param quad - velocity quadrature with abscissa and weights
   */
 
  
 DistFunctFields(const MeshList& meshes,const MacroFields& macroPr, const Quadrature<T>& quad):
  _meshes(meshes),
    _quad(quad) 
    {
      /**
       * integer N123
       * total number of velocity directions.
       */
      const int numFields = _quad.getDirCount();
      for(int n=0; n<numFields; n++)
	{
	  stringstream ss;
	  ss << n;
	  string fieldName = "dist_Field_" + ss.str(); 
	  dsf.push_back(new Field(fieldName));
	}
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& cells = mesh.getCells();
	  const int nCells = cells.getSelfCount();
	  double pi(3.14159);
	  
	  const TArray& density = dynamic_cast<const TArray&>(macroPr.density[cells]);
	  const TArray& temperature = dynamic_cast<const TArray&>(macroPr.temperature[cells]);
	  const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(macroPr.velocity[cells]);
	  
	  
	  const TArray& cx = dynamic_cast<const TArray&>(*_quad.cxPtr);
	  const TArray& cy = dynamic_cast<const TArray&>(*_quad.cyPtr);
	  const TArray& cz = dynamic_cast<const TArray&>(*_quad.czPtr);
	  
	  /*
	    TArray* distfunAPtr;   
	    distfunAPtr=new TArray(nCells);
	    TArray & distfunA= *distfunAPtr;  
	  */      
	  
	  for(int j=0;j<numFields;j++){
	    Field& fnd= *dsf[j]; 
	    
	    shared_ptr<TArray> fcPtr(new TArray(cells.getCount()));
	    
	    fnd.addArray(cells,fcPtr);
	    
	    //TArray& fc = dynamic_cast<TArray&>(fnd[cells]);
	    TArray& fc = *fcPtr;
	    for(int c=0; c<nCells;c++) {
		fc[c]=density[c]/pow((pi*temperature[c]),1.5)*
		  exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
			pow((cz[j]-v[c][2]),2.0))/temperature[c]);
	      } 
	    
	  }
	}
    }
  
 DistFunctFields(const MeshList& meshes, const Quadrature<T>& quad):
  _meshes(meshes),
    _quad(quad)
    {
      FILE * pFile;
      pFile=fopen("distfun.txt","w");
	/**
       * integer N123
       * total number of velocity directions.
       */
      const int numFields = _quad.getDirCount(); 
      for(int j=0; j<numFields; j++)
	{
	  stringstream ss;
	  ss << j;
	  string fieldName = "dist_Field_" + ss.str(); 
	  dsf.push_back(new Field(fieldName));
	}
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& cells = mesh.getCells();
	  const int nCells = cells.getSelfCount();
	  double pi(3.14159);
	  
	  const TArray& cx = dynamic_cast<const TArray&>(*_quad.cxPtr);
	  const TArray& cy = dynamic_cast<const TArray&>(*_quad.cyPtr);
	  const TArray& cz = dynamic_cast<const TArray&>(*_quad.czPtr);
	  const TArray& dcxyz = dynamic_cast<const TArray&>(*_quad.dcxyzPtr);
	  for(int j=0;j<numFields;j++){
	    Field& fnd= *dsf[j]; 
	    shared_ptr<TArray> fcPtr(new TArray(nCells));
	    
	    fnd.addArray(cells,fcPtr);
	    
	    TArray& fc = dynamic_cast<TArray&>(fnd[cells]);
	    //TArray& fc = *fcPtr;
	    
	    for(int c=0; c<nCells;c++){
	      fc[c]=1./pow(pi*1.0,1.5)*exp(-(pow((cx[j]-1.0),2.0)+pow((cy[j]-0.5),2.0)+					       pow((cz[j]-0.),2.0))/1.0);
	
	      } 
	    fprintf(pFile,"%12.6f %12.6f %12.6f %12.6f %E \n",cx[j],cy[j],cz[j],dcxyz[j],fc[0]);
	  }
	  fclose(pFile);
	}
    }

  /* 
 void init()
  { std::vector<Field*> dsf(new std::vector<Field*>(_quad.getDirCount()));
    for(int j=0;j<N123;j++){
      // Field& fnd;
      const int numMeshes = _meshes.size();
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  FlowVC<T> *vc(new FlowVC<T>());
	  vc->vcType = "flow";
	  _vcMap[mesh.getID()] = vc;
	  const FlowVC<T>& vc = *_vcMap[mesh.getID()];
	  
	  const StorageSite& cells = mesh.getCells();
	  const StorageSite& faces = mesh.getFaces();
	  shared_ptr<TArray> rhoCell(new TArray(cells.getCount()));
	  *rhoCell = vc["density"];
	  _macroFields.density.addArray(cells,rhoCell);
	}
    }
  }
  */
  void initializeMaxwellian(const MacroFields& macroPr, DistFunctFields<T>& dsfPtr)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();
	double pi(3.14159);
	
	const TArray& density = dynamic_cast<const TArray&>(macroPr.density[cells]);
	const TArray& temperature = dynamic_cast<const TArray&>(macroPr.temperature[cells]);
	const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(macroPr.velocity[cells]);
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quad.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quad.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quad.czPtr);
	const int numFields= _quad.getDirCount(); 
	
	for(int j=0;j< numFields;j++){
	  Field& fnd = *dsfPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=density[c]/pow((pi*temperature[c]),1.5)*
	      exp(-(pow((cx[j]-v[c][0]),2.0)+pow((cy[j]-v[c][1]),2.0)+
		    pow((cz[j]-v[c][2]),2.0))/temperature[c]);
	  } 
	  
	}
      }
  }
  
  
  
  
 private:
  const MeshList _meshes;
  const Quadrature<T> _quad;
  //FlowBCMap _bcMap;
  //FlowVCMap _vcMap;
};
 

/*
//Fgamma, Macroparameters, -ESBGKModel.h
void Fgamma_alphas(){
  for (int j=0;j<N123;j++){
    malpha_BGK[j][0]=1.0;
    malpha_BGK[j][1]=cx[j];
    malpha_BGK[j][2]=pow(cx[j],2.0)+pow(cy[j],2.0)+pow(cz[j],2.0);  
    malpha_BGK[j][3]=cy[j];
  }
  if(solver == 2){
    for (int j=0;j<N123;j++){
      malpha_ES[j][0]=1.0;
      malpha_ES[j][1]=cx[j];
      malpha_ES[j][2]=cy[j];
      malpha_ES[j][3]=pow(cx[j],2.0);
      malpha_ES[j][4]=pow(cy[j],2.0);
      malpha_ES[j][5]=pow(cz[j],2.0);
      malpha_ES[j][6]=cx[j]*cy[j];
    }}

}
*/

#endif
