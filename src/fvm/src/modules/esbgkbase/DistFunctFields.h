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
    _quadrature(quad) 
    {
      /**
       * integer N123
       * total number of velocity directions.
       */
      const int numFields = _quadrature.getDirCount();
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
	  
	  
	  const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	  const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	  const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	  
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
    _quadrature(quad)
    {
      FILE * pFile;
      pFile=fopen("distfun.txt","w");
	/**
       * integer N123
       * total number of velocity directions.
       */
      const int numFields = _quadrature.getDirCount(); 
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
	  const int nCells = cells.getCount();
	  double pi(3.14159);
	 
	  const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	  const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	  const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	  const TArray& dcxyz = dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
	  for(int j=0;j<numFields;j++){
	    Field& fnd= *dsf[j]; 
	    shared_ptr<TArray> fcPtr(new TArray(nCells));
	    
	    fnd.addArray(cells,fcPtr);
	    
	    TArray& fc = dynamic_cast<TArray&>(fnd[cells]);
	    //TArray& fc = *fcPtr;
	    
	    for(int c=0; c<nCells;c++){
	      fc[c]=1./pow(pi*1.0,1.5)*exp(-(pow((cx[j]-1.0),2.0)+pow((cy[j]-0.0),2.0)+					       pow((cz[j]-0.0),2.0))/1.0);
	
	      } 
	    fprintf(pFile,"%12.6f %12.6f %12.6f %12.6f %E \n",cx[j],cy[j],cz[j],dcxyz[j],fc[0]);
	  }
	  fclose(pFile);
	}
    }


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
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
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
  void weightedMaxwellian(DistFunctFields<T>& dsfPtr)
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int nCells = cells.getSelfCount();
	double pi(3.14159);
	
	const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
	const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
	const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
	const int numFields= _quadrature.getDirCount(); 
	
	for(int j=0;j< numFields;j++){
	  Field& fnd = *dsfPtr.dsf[j];
	  TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	  for(int c=0; c<nCells;c++){
	    f[c]=0.5*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-4.0),2.0)+pow((cy[j]-0.0),2.0)+pow((cz[j]-0.0),2.0))/1.0)
	      +0.5*1.0/pow((pi*1.0),1.5)*exp(-(pow((cx[j]-4.0),2.0)+pow((cy[j]-0.0),2.0)+pow((cz[j]-0.0),2.0))/1.0);
	  }
	}
      }
  }
  
  void OutputDistributionFunction(DistFunctFields<T>& dsfPtr) //, const char* filename)
  {
    FILE * pFile;
    pFile = fopen("outputf0.plt","w");  
    int N1=_quadrature.getNVCount();
    int N2=_quadrature.getNthetaCount();
    int N3=_quadrature.getNphiCount();
    fprintf(pFile,"%s \n", "VARIABLES= 'cx', 'cy', 'cz', 'f',");
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
	Field& fnd = *dsfPtr.dsf[j];
	TArray& f = dynamic_cast< TArray&>(fnd[cells]);
	fprintf(pFile,"%E %E %E %E\n",cx[j],cy[j],cz[j],f[0]);
      }
    }
  }

  
  
  
 private:
  const MeshList _meshes;
  const Quadrature<T> _quadrature;
  //KineticModelOptions<T> _options;
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
