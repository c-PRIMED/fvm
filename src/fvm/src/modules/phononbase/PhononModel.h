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
#include "PhononInterface.h"
#include "LinearSystem.h"
#include "NumType.h"
#include "PhononCollisionDiscretization.h"
#include "PhononConvectionDiscretization.h"
#include "GenericPhononBCS.h"
#include "Linearizer.h"
#include "ArrowHeadMatrix.h"
#include "COMETDiscretizer.h"
#include <math.h>

template<class T>
class PhononModel : public Model
{

 public:
  
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef shared_ptr<VectorT3Array> T3ptr;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef Array<T> Tarray;
  typedef shared_ptr<Tarray> Tarrptr;
  typedef map<int,PhononBC<T>*> PhononBCMap;
  typedef typename Tmode::Mode_ptr Mode_ptr;
  typedef typename Tmode::Reflection Reflection;
  typedef typename Tmode::Reflptr Reflptr;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef typename Tmode::Refl_Map Refl_Map;
  typedef Array<int> BCcellArray;
  typedef shared_ptr<BCcellArray> BCellPtr;
  typedef vector<BCellPtr> BCcellList;
  typedef Array<bool> BCfaceArray;
  typedef shared_ptr<BCfaceArray> BfacePtr;
  typedef vector<BfacePtr> BCfaceList;
  
 PhononModel(const MeshList& meshes,const GeomFields& geomFields,Tkspace& kspace,PhononMacro& macro):
  
  Model(meshes),
    _geomFields(geomFields),
    _kspace(kspace),
    _macro(macro),
    _BCells(),
    _BFaces()
      {	
      //setting boundary conditions
	const int numMeshes = _meshes.size();
  
      for (int n=0; n<numMeshes; n++) //mesh loop beg
	{
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& cells=mesh.getCells();
	  const StorageSite& faces=mesh.getFaces();
	  const int cellCount=cells.getSelfCount();
	  const int faceCount=faces.getCount();

	  BCellPtr BCptr(new BCcellArray(cellCount));
	  _BCells.push_back(BCptr);
	  BCptr->zero();

	  BfacePtr BFptr(new BCfaceArray(faceCount));
	  BFptr->zero();
	  _BFaces.push_back(BFptr);
	  
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
	    }
	}
    
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
	const StorageSite& cells=mesh.getCells();
	const int numcells=cells.getCount();
	
	//initialize lattice temperature
	shared_ptr<Tarray> TLcell(new Tarray(numcells));
	const T Tinit=_options["initialTemperature"];
	*TLcell=Tinit;
	_macro.temperature.addArray(cells,TLcell);

	T e0sum=0.;
	
	for (int k=0;k<numK;k++)  //kspace loop beg
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int numM=kv.getmodenum();
	    const T dk3=kv.getdk3();
	    
	    for (int m=0;m<numM;m++) //mode loop beg
	      {
		//initialize each phonon mode
		Tmode& mode=kv.getmode(m);
		Field& efield=mode.getfield();
		Field& e0field=mode.gete0field();
		Field& resfield=mode.getresid();
		Tarrptr e0var=shared_ptr<Tarray>(new Tarray(numcells));
		Tarrptr evar=shared_ptr<Tarray>(new Tarray(numcells));
		Tarrptr resid=shared_ptr<Tarray>(new Tarray(numcells));
		const T einit=mode.calce0(Tinit);
		e0sum+=einit*dk3;
		*evar=einit;
		*e0var=einit;
		*resid=0.;
		efield.addArray(cells,evar);
		e0field.addArray(cells,e0var);
		resfield.addArray(cells,resid);
	      }; //mode loop end
	  }; //kspace loop end

	e0sum=e0sum/_kspace.getDK3();
	shared_ptr<Tarray> e0cell(new Tarray(numcells));
	*e0cell=e0sum;
	_macro.e0.addArray(cells,e0cell);
	shared_ptr<Tarray> e0ResidCell(new Tarray(numcells));
	*e0ResidCell=0.;
	_macro.e0Residual.addArray(cells,e0ResidCell);

	// Compute reflections--later, should be moved inside above loop
	
	for (int k=0;k<numK;k++)  //kspace loop beg
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int numM=kv.getmodenum();
	    const T dk3=kv.getdk3();
	    
	    for (int m=0;m<numM;m++) //mode loop beg
	      {	
		Tmode& mode=kv.getmode(m);
		Refl_Map& rmap=mode.getreflmap();
		VectorT3 vg=mode.getv();
		T vmag=sqrt(pow(vg[0],2)+pow(vg[1],2)+pow(vg[2],2));
		VectorT3 si=vg/vmag;
		VectorT3 so;

		foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    const StorageSite& faces = fg.site;
		    const Field& AreaMagField=_geomFields.areaMag;
		    const Tarray& AreaMag=dynamic_cast<const Tarray&>(AreaMagField[faces]);
		    const Field& AreaDirField=_geomFields.area;
		    const VectorT3Array& AreaDir=
		      dynamic_cast<const VectorT3Array&>(AreaDirField[faces]);
		   	      
		    const VectorT3 n=AreaDir[0]/AreaMag[0];
		    const T sidotn=si[0]*n[0]+si[1]*n[1]+si[2]*n[2];

		    if (sidotn > T_Scalar(0.0))
		      {
			so=si-2.*(si[0]*n[0]+si[1]*n[1]+si[2]*n[2])*n;
			Refl_pair refls;
			_kspace.findSpecs(dk3,vmag,m,so,refls);
			rmap[fg.id]=refls;
		      }
		  }
	      }
	  }

	BCcellArray& BCArray=*(_BCells[n]);
	BCfaceArray& BCfArray=*(_BFaces[n]);
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if((_bcMap[fg.id]->bcType == "reflecting"))
	      {
		const CRConnectivity& BfaceCells=mesh.getFaceCells(fg.site);
		const int faceCount=fg.site.getCount();
		const int offSet=fg.site.getOffset();
		for(int i=0;i<faceCount;i++)
		  {
		    int cell1=BfaceCells(i,0);
		    BCArray[cell1]=1;
		  }
		for(int i=offSet;i<offSet+faceCount;i++)
		  BCfArray[i]=true;
	      } 
	  }	
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
	    const PhononBC<T>& bc = *_bcMap[fg.id];
	    
	    PhononBoundary<T> pbc(faces, mesh,_geomFields,_kspace,_options,fg.id);
	    
	    if (bc.bcType == "reflecting")
	      {	
     
		FloatValEvaluator<T>
		  bReflection(bc.getVal("specifiedReflection"),faces);

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
	foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const Mesh& otherMesh = faces.getMesh();
	    cout << " " << endl;
	    //cout << mesh << " " << otherMesh << endl;
	    //PhononInterface<T> pbc(faces, mesh,_geomFields,_kspace,_options,fg.id);
	  };
      };
  };
  
  void updateTL()
  {

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	const int numK=_kspace.getlength();
	Tarray& TL=dynamic_cast<Tarray&>(_macro.temperature[cells]);
	Tarrptr e_sumptr=shared_ptr<Tarray>(new Tarray(numcells));
	Tarray& e_sum=*(e_sumptr);
	e_sum=0.;
	
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
		    e_sum[c]+=e_val[c]*dk3/tau;

	      }
	  }
	
	for(int c=0;c<numcells;c++)
	  _kspace.NewtonSolve(TL[c],e_sum[c]);
      }
  }

  void COMETupdateTL()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	Tarray& TL=dynamic_cast<Tarray&>(_macro.temperature[cells]);
	Tarray& e0Array=dynamic_cast<Tarray&>(_macro.e0[cells]);
	
	for(int c=0;c<numcells;c++)
	  _kspace.NewtonSolve(TL[c],e0Array[c]);
      }
  }

  void updatee0()
  {
    
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	const int numK=_kspace.getlength();
	Tarray& TL=dynamic_cast<Tarray&>(_macro.temperature[cells]);
	
	for(int k=0;k<numK;k++)
	  {
	    Tkvol& kv=_kspace.getkvol(k);
	    const int modenum=kv.getmodenum();
	    
	    for(int m=0;m<modenum;m++)
	      {		
		Tmode& mode=kv.getmode(m);
		Field& e0field=mode.gete0field();
		Tarray& e0_val=dynamic_cast<Tarray&>(e0field[cells]);
		for(int c=0;c<numcells;c++)
		  e0_val[c]=mode.calce0(TL[c]);
	      }
	  }
	
      }
  }

 void updateHeatFlux()
 {
   const int numMeshes = _meshes.size();
   for (int n=0; n<numMeshes; n++)
     {
       VectorT3 zero_vec;
       const Mesh& mesh = *_meshes[n];
       const StorageSite& cells = mesh.getCells();
       const int numcells = cells.getCount();
       T3ptr heatFluxptr=T3ptr(new VectorT3Array(numcells));
       VectorT3Array& heatFlux=*heatFluxptr;
       zero_vec[0]=0.;zero_vec[1]=0.;zero_vec[2]=0.;
       heatFlux=zero_vec;
       
       const int numK=_kspace.getlength();
       for(int k=0;k<numK;k++)
	 {
	   Tkvol& kv=_kspace.getkvol(k);
	   T dk3=kv.getdk3();
	   const int modenum=kv.getmodenum();
	   for(int m=0;m<modenum;m++)
	     {
	       Tmode& mode=kv.getmode(m);
	       VectorT3 vg=mode.getv();
	       Field& efield=mode.getfield();
	       const Tarray& eval=dynamic_cast<const Tarray&>(efield[cells]);

	       for(int c=0;c<numcells;c++)
		 heatFlux[c]+=eval[c]*vg*dk3;
	     }
	 }
       _macro.heatFlux.addArray(cells,heatFluxptr);  
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

    int early=0;

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
	if(_niters % _options.showResidual == 0)
	  cout << _niters << ": " << *rNorm <<endl;

	_niters++;
	if ((*rNorm < _options.absTolerance)||(*normRatio < _options.relTolerance )) 
	  {
	    early=1;
	    cout << endl;
	    cout << "Final Residual at "<<_niters<<": "<< *rNorm <<endl;
	    break;}
	
	updateTL();	//update macroparameters
	updatee0();
	callBoundaryConditions();
      }

    updateHeatFlux();	  
  }

  void advanceCOMET(const int niters)
  {  
    
    T tol=_options.absTolerance;
    T aveResid=-1;
    T residChange=1;
    int exit=0;

    for(int n=0;n<niters;n++)
      {
	if(exit==1)
	  {
	    break;
	  }

	const int numMeshes=_meshes.size();
	for(int msh=0;msh<numMeshes;msh++)
	  {
	    const Mesh& mesh=*_meshes[msh];
	    const BCcellArray& BCArray=*(_BCells[msh]);
	    const BCfaceArray& BCfArray=*(_BFaces[msh]);
	    COMETDiscretizer<T> CDisc(mesh,_geomFields,_macro,
				      _kspace,_bcMap,BCArray,BCfArray);

	    CDisc.setfgFinder();
	    CDisc.findResid();
	    
	    if(aveResid==-1)
	      aveResid=CDisc.getAveResid();
	    else
	      {
		residChange=fabs(aveResid-CDisc.getAveResid()/aveResid);
		aveResid=CDisc.getAveResid();
	      }

	    if(n%_options.showResidual==0)
	      cout<<"[Iteration: "<<n<<" || Residual Change: "<<residChange<<
		" || Average Residual: "<<aveResid<<"]"<<endl;
	    
	    if(aveResid>tol)
	      {
		CDisc.COMETSolve();
		callBoundaryConditions();
	      }
	    else
	      {
		exit=1;
		break;
	      }
	  }
      }
    COMETupdateTL();
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
	   cout<<TL[i]<<" "<<i<<endl;
       }

  }

  T HeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    T r(0.);
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
        {
	  const StorageSite& faces = fg.site;
	  const int nFaces = faces.getCount();
	  const StorageSite& cells = mesh.getCells();
	  const CRConnectivity& faceCells=mesh.getFaceCells(faces);
	  const Field& areaField=_geomFields.area;
	  const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]);
	    
	  for(int k=0;k<_kspace.getlength();k++)
	    {
	      Tkvol& kv=_kspace.getkvol(k);
	      int modenum=kv.getmodenum();
	      for(int m=0;m<modenum;m++)
		{
		  VectorT3 vg=kv.getmode(m).getv();
		  T dk3=kv.getdk3();
		  Field& efield=kv.getmode(m).getfield();
		  const Tarray& eval=dynamic_cast<const Tarray&>(efield[cells]);
		  for(int f=0; f<nFaces; f++)
		    {
		      const VectorT3 An=faceArea[f];
		      const int c1=faceCells(f,1);
		      const T vgdotAn=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
		      r += eval[c1]*vgdotAn*dk3;
		    }
		  found=true;
		}
	    }
        }
    }
    if (!found)
      throw CException("getHeatFluxIntegral: invalid faceGroupID");
    return r;
  }

  private:

    const GeomFields& _geomFields;
    Tkspace& _kspace;       //kspace
    PhononMacro& _macro;
    PhononModelOptions<T> _options;
    PhononBCMap _bcMap;
    MFRPtr _initialnorm;
    int _niters;
    BCcellList _BCells;
    BCfaceList _BFaces;
};

#endif
