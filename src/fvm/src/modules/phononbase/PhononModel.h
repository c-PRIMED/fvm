#ifndef _PHONONMODEL_H_
#define _PHONONMODEL_H_

#include "Model.h"
#include "Array.h"
#include "Vector.h"
#include <vector>
#include "Mesh.h"
#include "Kspace.h"
#include "kvol.h"
#include "pmode.h"
#include "PhononMacro.h"
#include "PhononBC.h"
#include <map>
#include "PhononBoundary.h"
#include "LinearSystem.h"

template<class T>
class PhononModel : public Model
{

 public:
  
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef Array<T> Tarray;
  typedef map<int,PhononBC<T>*> PhononBCMap;

 PhononModel(const MeshList& meshes,const GeomFields& geomFields,Tkspace& kspace,PhononMacro& macro):
  
  Model(meshes),
    _geomFields(geomFields),
    _kspace(kspace),
    _macro(macro)
      {	
      //setting boundary conditions
	const int numMeshes = _meshes.size();
  
      for (int n=0; n<numMeshes; n++) //mesh loop beg
	{
	  const Mesh& mesh = *_meshes[n];
	  
	  foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())  //facegroup loop beg
	    {
	      const FaceGroup& fg = *fgPtr;
	      if (_bcMap.find(fg.id) == _bcMap.end())
		{
		  PhononBC<T> *bc(new PhononBC<T>());
		  
		  _bcMap[fg.id] = bc;
		  if((fg.groupType == "wall"))
		    {
		      bc->bcType = "reflecting";
		    }
		  else if (fg.groupType == "velocity-inlet")
		    {
		      bc->bcType = "temperature";
		    }
		  else
		    throw CException("PhononModel: unknown face group type "
				     + fg.groupType);

		}
	    } //facegroup loop end
	} //mesh loop end
 
    };
  
  void init()
  {
    const int numMeshes=_meshes.size();

    for (int n=0;n<numMeshes;n++)  //mesh loop beg
      {
	const Mesh& mesh=*_meshes[n];
	const int numK=_kspace.getlength();
	const T DK3=_kspace.getDK3();     //total kspace volume
	const StorageSite& cells=mesh.getCells();
	const T Tref=_options["Tref"];
	const int numcells=cells.getSelfCount();
	
	//initialize lattice temperature
	shared_ptr<Tarray> TLcell(new Tarray(numcells));
	const T Tinit=_options["initialTemperature"];
	*TLcell=Tinit;
	_macro.temperature.addArray(cells,TLcell);

	for (int k=0;k<numK;k++)  //kspace loop beg
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int numM=kv.getmodenum();	    

	    for (int m=0;m<numM;m++) //mode loop beg
	      {
		//initialize each phonon mode
		Tmode& mode=kv.getmode(m);
		T cp=mode.getcp();
		Field& efield=mode.getfield();
		shared_ptr<Tarray> evar(new Tarray(numcells));
	    
		const T einit=cp*(Tinit-Tref)/DK3;
		*evar=einit;
		efield.addArray(cells,evar);

	      }; //mode loop end
	  }; //kspace loop end
      }; //mesh loop end
  };

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
	    const PhononBC<T>& bc = *_bcMap[fg.id];
	    
	    PhononBoundary<T> pbc(faces, mesh,_geomFields,_kspace,_options);

	    if (bc.bcType == "reflecting")
	      {	
     
		FloatValEvaluator<T>
		  bReflection(bc.getVal("specifiedreflection"),faces);

		pbc.applyReflectingWall(bReflection);
	      }
	    else if(bc.bcType=="temperature")
	      {
	      
		FloatValEvaluator<T>
		  bTemperature(bc.getVal("specifiedTemperature"),faces);
	  
		pbc.applyTemperatureWall(bTemperature);
	      } 
	  };
      };
  };

   void initPhononModelLinearization(LinearSystem& ls, Tmode& mode)
  {
   const int numMeshes = _meshes.size();
   for (int n=0; n<numMeshes; n++)
     {
       const Mesh& mesh = *_meshes[n];
       const StorageSite& cells = mesh.getCells();
       
       Field& fnd = mode.getfield();
       //const TArray& f = dynamic_cast<const TArray&>(fnd[cells]);
       
       MultiField::ArrayIndex vIndex(&fnd,&cells); //dsf is in 1 direction
       
       ls.getX().addArray(vIndex,fnd.getArrayPtr(cells));
       
       const CRConnectivity& cellCells = mesh.getCellCells();
       
       shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells)); //diagonal off diagonal, phi DiagTensorT3,T,VectorT3 for mometum in fvmbase 
       
       ls.getMatrix().addMatrix(vIndex,vIndex,m);
     }
   //mesh.getBoudaryfaces and interfacefaces?
  }; 

  /*

   void linearizeKineticModel(LinearSystem& ls, Tmode& mode)
   {
     // _velocityGradientModel.compute();
     DiscrList discretizations;
    
     Field& fnd = *mode.getfieldptr();
     T cp=mode.getcp();
     T tau=mode.gettau();
     
     const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
     const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
     const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
     
     if(_options.fgamma==2){ 
       Field& feqES = *_dsfEqPtrES.dsf[direction];
       shared_ptr<Discretization>
	 sdEQ(new CollisionTermDiscretization<T,T,T>
	      (_meshes, _geomFields, 
	       fnd,_macro.e0,tau); 
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
	     
	     if (( bc.bcType == "CopyBC")) 
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
  */   

  private:

    const GeomFields& _geomFields;
    Tkspace& _kspace;       //kspace
    PhononMacro& _macro;
    PhononModelOptions<T> _options;
    PhononBCMap _bcMap;

};

#endif
