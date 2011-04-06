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
#include "NumType.h"
#include "PhononCollisionDiscretization.h"
#include "PhononConvectionDiscretization.h"
#include "GenericPhononBCS.h"
#include "Linearizer.h"

template<class T>
class PhononModel : public Model
{

 public:
  
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef Array<T> Tarray;
  typedef shared_ptr<Tarray> Tarrptr;
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
      
      init();
      cout<<"Model Initialized"<<endl;
    }

  PhononModelOptions<T>& getOptions() {return _options;}
  PhononBCMap& getBCs() {return _bcMap;}
  
  void init()
  {

    _initialnorm = MFRPtr();
    _niters=0;
    const int numMeshes=_meshes.size();
    
    for (int n=0;n<numMeshes;n++)  //mesh loop beg
      {
	const Mesh& mesh=*_meshes[n];
	const int numK=_kspace.getlength();
	const T DK3=_kspace.getDK3();     //total kspace volume
	const StorageSite& cells=mesh.getCells();
	const T Tref=_options["Tref"];
	const int numcells=cells.getCount();
	
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
		Field& e0field=mode.gete0field();		
		Tarrptr evar=shared_ptr<Tarray>(new Tarray(numcells));
		Tarrptr e0var=shared_ptr<Tarray>(new Tarray(numcells));		
		const T einit=cp*(Tinit-Tref)/DK3;
		*evar=einit;
		*e0var=einit;
		efield.addArray(cells,evar);
		e0field.addArray(cells,e0var);

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
	    else
	      {
		cout<<"Couldn't find boundary condition"<<endl;
		break;
	      }
	  };
      };
  };
  
  void updateTL()
  {

    T Tref=_options["Tref"];

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	const int numK=_kspace.getlength();
	const T Covertau=_kspace.getCovertau();
	Tarray& TL=dynamic_cast<Tarray&>(_macro.temperature[cells]);
	TL=0.;  // initailize to zero for summing
	
	for(int k=0;k<numK;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    T dk3=kv.getdk3();
	    
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		T tau=mode.gettau();
		Field& efield=mode.getfield();
		Tarray& e_val=dynamic_cast<Tarray&>(efield[cells]);
		
		for(int c=0;c<numcells;c++)
		    TL[c]+=e_val[c]*dk3/tau;

	      }
	  }
	
	for(int c=0;c<numcells;c++)
	    TL[c]=TL[c]/Covertau+Tref;

      }
  }

void updatee0()
  {

    T Tref=_options["Tref"];

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	const int numK=_kspace.getlength();
	T DK3=_kspace.getDK3();
	Tarray& TL=dynamic_cast<Tarray&>(_macro.temperature[cells]);
	
	for(int k=0;k<numK;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    
	    for(int m=0;m<modenum;m++)
	      {
		Tmode& mode=kv.getmode(m);
		T cp=mode.getcp();
		Field& e0field=mode.gete0field();
		Tarray& e0_val=dynamic_cast<Tarray&>(e0field[cells]);
		
		for(int c=0;c<numcells;c++)
		  e0_val[c]=cp*(TL[c]-Tref)/DK3;

	      }
	  }
	
      }
  }

   void initPhononModelLinearization(LinearSystem& ls, Tmode& mode)
  {
   const int numMeshes = _meshes.size();
   for (int n=0; n<numMeshes; n++)
     {
       const Mesh& mesh = *_meshes[n];
       const StorageSite& cells = mesh.getCells();
       
       Field& fnd = mode.getfield();
       
       MultiField::ArrayIndex vIndex(&fnd,&cells); 
       
       ls.getX().addArray(vIndex,fnd.getArrayPtr(cells));
       
       const CRConnectivity& cellCells = mesh.getCellCells();
       
       shared_ptr<Matrix> m(new CRMatrix<T,T,T>(cellCells));
       
       ls.getMatrix().addMatrix(vIndex,vIndex,m);
     }
   
  }; 
   
   void linearizePhononModel(LinearSystem& ls, Tmode& mode)
   {
     // _velocityGradientModel.compute();
     DiscrList discretizations;
    
     Field& fnd = mode.getfield();  //field for e"
     Field& e0fld=mode.gete0field(); //field for e0
     T tau=mode.gettau();
     VectorT3 vg=mode.getv();
    
     shared_ptr<Discretization>
       colldisc(new PhononCollisionDiscretization<T,T,T>(_meshes, _geomFields,fnd,e0fld,tau)); 
     
     discretizations.push_back(colldisc);
     
     shared_ptr<Discretization> 
       convdisc(new PhononConvectionDiscretization<T,T,T> (_meshes,_geomFields,fnd,vg));

     discretizations.push_back(convdisc);

     /*
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
	 }*/

     
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
	     const PhononBC<T>& bc = *_bcMap[fg.id];
	     Tarray& dsf = dynamic_cast< Tarray&>(fnd[cells]);
	     BaseGenericPhononBCS<T,T,T> gpbc(faces, mesh, _geomFields,
					       fnd,
					       ls.getMatrix(),
					       ls.getX(),
					       ls.getB());
	     
	     
	     const CRConnectivity& faceCells = mesh.getFaceCells(faces);			 
	     const Field& areaMagField = _geomFields.areaMag;
	     const Tarray& faceAreaMag = dynamic_cast<const Tarray &>(areaMagField[faces]);
	     const Field& areaField = _geomFields.area;
	     const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]); 
	     
	     if (( bc.bcType == "CopyBC")) 
	       {
		 for(int f=0; f< nFaces; f++)
		   {
		     gpbc.applyExtrapolationBC(f);
		   }
	       }     
	     else{
	       for(int f=0; f< nFaces; f++)
		 {
		   const VectorT3 en = faceArea[f]/faceAreaMag[f];
		   const T v_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
		   
		   if(v_dot_en < T_Scalar(0.0))
		     {		     
		       //incoming direction - dirchlet bc
		       const int c1= faceCells(f,1);    // boundary cell
		       T bvalue =dsf[c1];
		       gpbc.applyDirichletBC(f,bvalue);
		     }
		   else
		     {
		     //outgoing direction - extrapolation bc
		     gpbc.applyExtrapolationBC(f);
		     }    
		 }
	     }   
	   }
       }

   }
  
  void advance(const int niter)
  {
    for(int n=0; n<niter; n++)  
      {
	const int klength =_kspace.getlength();
	MFRPtr rNorm;
    
	/*
   	  _macroFields.velocity.syncLocal();
 	  _macroFields.temperature.syncLocal();   ???????
  	  _macroFields.density .syncLocal();
	*/	

	for(int k=0; k<klength;k++)
	  {
	    
	    Tkvol& kv=_kspace.getkvol(k);
	    const int mlength=kv.getmodenum();
	    
	    for(int m=0;m<mlength;m++)
	      {
		
		Tmode& mode=kv.getmode(m);
		
		LinearSystem ls;
		initPhononModelLinearization(ls,mode);
		ls.initAssembly();
		linearizePhononModel(ls,mode);
		ls.initSolve();
		
		//const T* kNorm=_options.getKineticLinearSolver().solve(ls);
		MFRPtr kNorm(_options.getPhononLinearSolver().solve(ls));
		
		
		if (!rNorm)
		  rNorm = kNorm;
		else
		  {
		    // find the array for the 0the direction residual and
		    // add the current residual to it
		    Tkvol& kvol0=_kspace.getkvol(0);
		    Tmode& mode0=kvol0.getmode(0);
		    Field& efield0=mode0.getfield();
		    Field& efield=mode.getfield();
		    
		    ArrayBase& rArray0 = (*rNorm)[efield0];
		    ArrayBase& rArrayd = (*kNorm)[efield];
		    rArray0 += rArrayd;
		  }
		
		ls.postSolve();
		ls.updateSolution();
	    	
		_options.getPhononLinearSolver().cleanup();
		
	      }
	  }
	
	
	if (!_initialnorm) _initialnorm = rNorm;

	if (_niters < 5)
	   _initialnorm->setMax(*rNorm);
            
 
	MFRPtr normRatio((*rNorm)/(*_initialnorm));	
	//	MFRPtr vnormRatio((*vNorm)/(*_initialKmodelvNorm));
	//if ( MPI::COMM_WORLD.Get_rank() == 0 )
	{cout << _niters << ": " << *rNorm <<endl; }

	_niters++;
	if ((*rNorm < _options.absTolerance)||(*normRatio < _options.relTolerance )) 
	  {//&& ((*vNorm < _options.absoluteTolerance)||(*vnormRatio < _options.relativeTolerance )))
	    break;}
	
	callBoundaryConditions();
	updateTL();	//update macroparameters
	updatee0();

      }
  }

  void printTemp()
  {
     const int numMeshes = _meshes.size();
     for (int n=0; n<numMeshes; n++)
       {
	 const Mesh& mesh = *_meshes[n];
	 const StorageSite& cells = mesh.getCells();
	 const int numcells=cells.getCount();
	 Tarray& TL=dynamic_cast<Tarray&>(_macro.temperature[cells]);
	 for(int i=0;i<numcells;i++)
	   cout<<TL[i]<<endl;

       }

  }    

  private:

    const GeomFields& _geomFields;
    Tkspace& _kspace;       //kspace
    PhononMacro& _macro;
    PhononModelOptions<T> _options;
    PhononBCMap _bcMap;
    MFRPtr _initialnorm;
    int _niters;
    

};

#endif
