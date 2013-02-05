// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _COMETMODEL_H_
#define _COMETMODEL_H_

#include "Model.h"
#include "Array.h"
#include "Mesh.h"
#include "Kspace.h"
#include "kvol.h"
#include "pmode.h"
#include "COMETDiscretizer.h"
#include "NumType.h"
#include "COMETBC.h"
#include "PhononMacro.h"
#include "COMETBoundary.h"
#include "COMETInterface.h"

template<class T>
class COMETModel : public Model
{

 public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef shared_ptr<VectorT3Array> VT3Ptr;
  typedef Kspace<T> Tkspace;
  typedef Tkspace* TkspPtr;
  typedef vector<Tkspace*> TkspList;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef Array<bool> BArray;
  typedef Array<int> IntArray;
  typedef shared_ptr<IntArray> IntArrPtr;
  typedef Array<T> TArray;
  typedef shared_ptr<TArray> TArrptr;
  typedef map<int,COMETBC<T>*> COMETBCMap;
  typedef vector<COMETIC<T>*> IClist;
  typedef map<int,int> MeshKspaceMap;
  typedef COMETModel<T> TCOMET;
  typedef shared_ptr<TCOMET> TCOMETPtr;
  typedef COMETModelOptions<T> TCModOpts;
  typedef typename Tmode::Mode_ptr Mode_ptr;
  typedef typename Tmode::Reflection Reflection;
  typedef typename Tmode::Reflptr Reflptr;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef typename Tmode::Refl_Map Refl_Map;
  typedef shared_ptr<MeshList> MshLstPtr;
  typedef shared_ptr<Mesh> MeshPtr;
  typedef shared_ptr<GeomFields> GeoFldsPtr;
  typedef shared_ptr<Tkspace> KspPtr;
  typedef shared_ptr<PhononMacro> PMacroPtr;
  typedef shared_ptr<StorageSite> SSPtr;
  typedef shared_ptr<CRConnectivity> CRPtr;
  typedef Array<int> BCfaceArray;
  typedef shared_ptr<BCfaceArray> BfacePtr;
  typedef vector<BfacePtr> BCfaceList;
  typedef Array<int> BCcellArray;
  typedef shared_ptr<BCcellArray> BCellPtr;
  typedef vector<BCellPtr> BCcellList;
  typedef vector<IntArrPtr> MeshICmap;
  typedef KSConnectivity<T> TKConnectivity;
  typedef TKConnectivity* TKCptr;
  typedef vector<TKCptr> TKClist;
  typedef map<int, TKClist> FgTKClistMap;
  typedef vector<shared_ptr<Field> > FieldVector;
  typedef pair<const StorageSite*, const StorageSite*> SSPair;
  typedef map<SSPair,shared_ptr<Array<int> > > ScatGathMaps;
  typedef map<const StorageSite*, StorageSite*> SiteMap;   //fine to coarse
  
 COMETModel(const MeshList& meshes, const int level, GeomFields& geomFields,
	    TkspList& kspaces, PhononMacro& macro):
  Model(meshes),
    _level(level),
    _geomFields(geomFields),
    _kspaces(kspaces),
    _macro(macro),
    _finestLevel(NULL),
    _coarserLevel(NULL),
    _finerLevel(NULL),
    _residual(1.0),
    _MeshToIC(),
    _rank(-1)
      {
	const int numMeshes = _meshes.size();

#ifdef FVM_PARALLEL
	_rank=MPI::COMM_WORLD.Get_rank();
#endif
	
	for (int n=0; n<numMeshes; n++)
	  {
	    Mesh& mesh = *_meshes[n];
	    const StorageSite& faces=mesh.getFaces();
	    const StorageSite& cells=mesh.getCells();
	    const int faceCount=faces.getCount();
	    const int cellCount=cells.getCount();
	    _MeshKspaceMap[n]=-1; //mesh map not set
	    
	    BfacePtr BFptr(new BCfaceArray(faceCount));
	    BFptr->zero();
	    _BFaces.push_back(BFptr);
	    
	    BCellPtr BCptr(new BCcellArray(cellCount));
	    _BCells.push_back(BCptr);
	    BCptr->zero();
	    
	    if(_level==0)
	      {
		foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    if (_bcMap.find(fg.id) == _bcMap.end() && fg.id>0)
		      {
			COMETBC<T> *bc(new COMETBC<T>()); 
			_bcMap[fg.id] = bc;
		      }
		  }
	      }
	  }
      }
  
  TCModOpts& getOptions() {return _options;}
  COMETBCMap& getBCs() {return _bcMap;}
  MeshKspaceMap& getMKMap() {return _MeshKspaceMap;}

  void init()
  {
    const int numMeshes=_meshes.size();
    const T Tinit=_options["initialTemperature"];

     for (int n=0;n<numMeshes;n++)
       {
	 Mesh& mesh=*_meshes[n];
	 if(_MeshKspaceMap[n]==-1)
	   throw CException("Have not set the Kspace for this Mesh!!");
	 Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	 const int kcount=kspace.gettotmodes();
	 BCfaceArray& BCfArray=*(_BFaces[n]);
	 BCcellArray& BCArray=*(_BCells[n]);
	 const int numK=kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 VectorT3 lamTemp;

	 // set micro parameters
	 cout<<"Allocating Arrays..."<<endl;
	 shared_ptr<TArray> eArray(new TArray(numcells*kcount));
	 shared_ptr<TArray> e0Array(new TArray(numcells*kcount));
	 shared_ptr<TArray> ResidArray(new TArray(numcells*kcount));
	 shared_ptr<TArray> tauArray(new TArray(numcells*kcount));
	 kspace.seteArray(eArray);
	 kspace.sete0Array(e0Array);
	 kspace.setResArray(ResidArray);
	 kspace.setTauArray(tauArray);
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 shared_ptr<TArray> deltaTcell(new TArray(numcells));
	 shared_ptr<VectorT3Array> lamArray(new VectorT3Array(numcells));
	 shared_ptr<VectorT3Array> q(new VectorT3Array(numcells));
	 lamTemp.zero();
	 q->zero();
	 *lamArray=lamTemp;
	 *deltaTcell=0.;
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 _macro.deltaT.addArray(cells,deltaTcell);
	 _macro.lam.addArray(cells,lamArray);
	 _macro.heatFlux.addArray(cells,q);

	 shared_ptr<IntArray> f2c(new IntArray(numcells));
	 *f2c=-1;
	 _geomFields.fineToCoarse.addArray(cells, f2c);

	 int modeCount=kspace.getkvol(0).getmodenum();
	 FieldVector* FieldVecPtr=new FieldVector();

	 for(int mode=0;mode<modeCount;mode++)
	   {
	     shared_ptr<Field> modeField(new Field("mode"));
	     shared_ptr<TArray> modeTemp(new TArray(numcells));
	     *modeTemp=Tinit;
	     modeField->addArray(cells,modeTemp);
	     FieldVecPtr->push_back(modeField);
	   }

	 _macro.BranchTemperatures[n]=FieldVecPtr;

	 cout<<"Arrays Allocated...Initializing Values..."<<endl;
	 
	 for(int c=0;c<numcells;c++)
	   {
	     int cellIndex=kspace.getGlobalIndex(c,0);
	     for (int k=0;k<numK;k++)
	       {
		 Tkvol& kv=kspace.getkvol(k);
		 const int numM=kv.getmodenum();
		 const T dk3=kv.getdk3();

		 for (int m=0;m<numM;m++)
		   {
		     Tmode& mode=kv.getmode(m);
		     T tau=mode.gettau();
		     const T einit=mode.calce0(Tinit);
		     /*
		     if(m==0 && k==0)
		       (*eArray)[cellIndex]=1.1*einit;
		       else*/

		     (*eArray)[cellIndex]=einit;
		     (*e0Array)[cellIndex]=einit;
		     (*ResidArray)[cellIndex]=0.;
		     (*tauArray)[cellIndex]=tau;
		     cellIndex++;
		   }
	       }
	     kspace.updateTau(c,Tinit);
	   }
	 
	 
	 shared_ptr<TArray> TlResidCell(new TArray(numcells));
	 *TlResidCell=0.;
	 _macro.TlResidual.addArray(cells,TlResidCell);

	 cout<<"Values Initialized...Setting Facegroups..."<<endl;
	 
	 //setting facegroups

	 foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     const StorageSite& faces = fg.site;
	     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);

	     const int faceCount=faces.getCount();
	     const int offSet=faces.getOffset();

	     for(int i=offSet;i<offSet+faceCount;i++)
	       BCfArray[i]=-1;  //implicit boundary
	   }
	 

	 foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     if(fg.id>0)
	       {
		 if(_bcMap[fg.id]->bcType == "reflecting")
		   {
		     const StorageSite& faces = fg.site;
		     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
		     const int faceCount=faces.getCount();
		     const int offSet=faces.getOffset();
		     const bool Imp=(*(_bcMap[fg.id]))["FullyImplicit"];
		     const T ref=(*(_bcMap[fg.id]))["specifiedReflection"];
		     const Field& AreaMagField=_geomFields.areaMag;
		     const TArray& AreaMag=
		       dynamic_cast<const TArray&>(AreaMagField[faces]);
		     const Field& AreaDirField=_geomFields.area;
		     const VectorT3Array& AreaDir=
		       dynamic_cast<const VectorT3Array&>(AreaDirField[faces]);

		     const VectorT3 norm=AreaDir[0]/AreaMag[0];

		     if(ref==1.)
		       {
			 for (int k=0;k<numK;k++)
			   {
			     Tkvol& kv=kspace.getkvol(k);
			     const int numM=kv.getmodenum();
			     const T dk3=kv.getdk3();
			     for (int m=0;m<numM;m++)
			       {

				 Tmode& mode=kv.getmode(m);
				 const VectorT3 vg=mode.getv();
				 const T vmag=sqrt(pow(vg[0],2)+pow(vg[1],2)+pow(vg[2],2));
				 VectorT3 si=vg/vmag;
				 VectorT3 so;
				 const T sidotn=si[0]*norm[0]+si[1]*norm[1]+si[2]*norm[2];
				 Refl_Map& rmap=mode.getreflmap();

				 if (sidotn > T_Scalar(0.0))
				   {
				     so=si-2.*(si[0]*norm[0]+si[1]*norm[1]+si[2]*norm[2])*norm;
				     T soMag=sqrt(pow(so[0],2)+pow(so[1],2)+pow(so[2],2));
				     so/=soMag;
				     so*=vmag;
				     Refl_pair refls;
				     Refl_pair reflsFrom;
				     kspace.findSpecs(norm,m,k,refls);
				     rmap[fg.id]=refls;
				     const int k1=refls.first.second;
				     Tmode& mode2=kspace.getkvol(k1).getmode(m);
				     Refl_Map& rmap2=mode2.getreflmap();
				     reflsFrom.first.second=-1;
				     reflsFrom.second.second=k;
				     rmap2[fg.id]=reflsFrom;
				   }
			       }
			   }
		       }
		     
		     if(Imp)
		       {
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=2;  //implicit boundary
		       }
		     else
		       {
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=3;  //explicit boundary
		       }
		     
		     if(Imp)
		       {
			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     if(BCArray[cell1]==0)
			       BCArray[cell1]=1;    //implicit boundary only
			     else if(BCArray[cell1]==2)
			       BCArray[cell1]=3;  //mix implicit/explicit boundary
			   }
		       }
		     else
		       {
			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     if(BCArray[cell1]==0)
			       BCArray[cell1]=2;  //explicit boundary only
			     else if (BCArray[cell1]==1)
			       BCArray[cell1]=3;  //mix implicit/explicit boundary
			   }
		       }
		   }
		 else if(_bcMap[fg.id]->bcType == "temperature")
		   {
		     const StorageSite& faces = fg.site;
		     const int faceCount=faces.getCount();
		     const int offSet=faces.getOffset();
		     
		     for(int i=offSet;i<offSet+faceCount;i++)
		       BCfArray[i]=1;
		   }
		 else if(_bcMap[fg.id]->bcType == "Interface")
		   {		 
		     StorageSite& faces0 = fg.site;
		     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces0);
		     const int faceCount=faces0.getCount();
		     const int offSet=faces0.getOffset();
		     
		     for(int i=offSet;i<offSet+faceCount;i++)
		       {
			 BCfArray[i]=4;  //Interface boundary
			 int cell0=BfaceCells(i-offSet,0);
			 int cell1=BfaceCells(i-offSet,1);
			 BCArray[cell0]=2;  //always treated explicitly
			 BCArray[cell1]=4; //interface ghost cells need to be labeled
		       }
		     
		     bool doneAlready=false;
		     foreach(const COMETIC<T>* icPtr, _IClist)
		       {
			 if(icPtr->FgID1==fg.id && icPtr->MeshID1==mesh.getID())
			   {
			     doneAlready=true;
			     break;
			   }
		       }
		     
		     bool foundMatch=false;
		     if(!doneAlready)
		       {
			 StorageSite* faces1Ptr=NULL;
			 int otherFgID,otherMid;
			 for(int otherMeshID=n+1;otherMeshID<numMeshes;otherMeshID++)
			   {
			     const Mesh& otherMesh=*_meshes[otherMeshID];
			     foreach(const FaceGroupPtr otherfgPtr, otherMesh.getBoundaryFaceGroups())
			       {
				 foundMatch=mesh.COMETfindCommonFaces(faces0, otherfgPtr->site, _geomFields);
				 if(foundMatch)
				   {
				     otherFgID=otherfgPtr->id;
				     otherMid=otherMeshID;
				     faces1Ptr=&(otherfgPtr->site);
				     break;
				   }
			       }
			     if(foundMatch)
			       break;
			   }
			 
			 if(foundMatch)
			   {
			     StorageSite& faces1=*faces1Ptr;
			     const Mesh& mesh1=*_meshes[otherMid];
			     COMETIC<T>* icPtr(new COMETIC<T>(n,fg.id,mesh1.getID(),
							      otherFgID, faces1.getCount()));
			     
			     icPtr->InterfaceModel=_bcMap[fg.id]->InterfaceModel;
			     
			     setLocalScatterMaps(mesh, faces0, mesh1, faces1);
			     
			     _IClist.push_back(icPtr);
			   }
			 else if(!doneAlready && !foundMatch)
			   {
			     cout<<"Face Group: "<<fg.id<<" MeshID: "<<mesh.getID()<<endl;
			     throw CException("Could not find a matching face group!");
			   }
		       }
		   }// end if interface 
	       }
	   }//end facegroup loop

	 cout<<"Facegroups Set...Mesh "<<n<<" Complete."<<endl;
	 
       }//end meshes loop

     cout<<"Mesh Loop Complete...Creating Interfaces (if any)..."<<endl;

     //Make map from mesh to IC
     IntArray ICcount(numMeshes);
     ICcount.zero();
     
     foreach(const COMETIC<T>* icPtr, _IClist)
       {
	 ICcount[icPtr->MeshID0]++;
	 ICcount[icPtr->MeshID1]++;
       }

     if(!_IClist.empty())
       {
	 for(int i=0;i<numMeshes;i++)
	   {
	     IntArrPtr MeshToIC(new IntArray(ICcount[i]));
	     _MeshToIC.push_back(MeshToIC);
	   }
       }
     
     ICcount.zero();
     const int listSize=_IClist.size();
     for(int ic=0;ic<listSize;ic++)
       {
	 const COMETIC<T>* icPtr=_IClist[ic];
	 (*(_MeshToIC[icPtr->MeshID0]))[ICcount[icPtr->MeshID0]]=ic;
	 (*(_MeshToIC[icPtr->MeshID1]))[ICcount[icPtr->MeshID1]]=ic;
	 ICcount[icPtr->MeshID0]++;
	 ICcount[icPtr->MeshID1]++;
       }

     COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);

     for(int ic=0;ic<listSize;ic++)
       {
	 COMETIC<T>* icPtr=_IClist[ic];
	 const int mid0=icPtr->MeshID0;
	 const int mid1=icPtr->MeshID1;
	 if(icPtr->InterfaceModel=="DMM")
	   ComInt.makeDMMcoeffs(*icPtr);
	 else if(icPtr->InterfaceModel=="NoInterface")
	   ComInt.makeNoInterfaceCoeffs(*icPtr);
	 ComInt.updateOtherGhost(*icPtr,mid0,false);
	 ComInt.updateOtherGhost(*icPtr,mid1,false);
       }

     cout<<"Interfaces Complete..."<<endl;
     
     initializeTemperatureBoundaries();
     _residual=updateResid(false);

     if(_options.DomainStats=="Loud")
       {
	 cout<<"Creating Coarse Levels on rank "<<_rank<<"..."<<endl;
	 cout<<"Level: 0, rank: "<<_rank<<endl<<endl;
	 calcDomainStats();
	 cout<<endl;
       }

     MakeCoarseModel(this);
     setFinestLevel(this);
     if(_options.DomainStats=="Loud")
       cout<<"Coarse Levels Completed on rank "<<_rank<<"."<<endl;
  }

  void initFromOld()
  {
    const int numMeshes=_meshes.size();

     for (int n=0;n<numMeshes;n++)
       {
	 Mesh& mesh=*_meshes[n];

	 BCfaceArray& BCfArray=*(_BFaces[n]);
	 BCcellArray& BCArray=*(_BCells[n]);
	 //const StorageSite& cells=mesh.getCells();
	 //const int numcells=cells.getCount();

	 //setting facegroups
	 foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     if(fg.id>0)
	       {
		 if(_bcMap[fg.id]->bcType == "reflecting")
		   {
		     const StorageSite& faces = fg.site;
		     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
		     const int faceCount=faces.getCount();
		     const int offSet=faces.getOffset();
		     const bool Imp=(*(_bcMap[fg.id]))["FullyImplicit"];
		     
		     if(Imp)
		       {
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=2;  //implicit boundary
		       }
		     else
		       {
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=3;  //explicit boundary
		       }
		     
		     if(Imp)
		       {
			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     if(BCArray[cell1]==0)
			       BCArray[cell1]=1;    //implicit boundary only
			     else if(BCArray[cell1]==2)
			       BCArray[cell1]=3;  //mix implicit/explicit boundary
			   }
		       }
		     else
		       {
			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     if(BCArray[cell1]==0)
			       BCArray[cell1]=2;  //explicit boundary only
			     else if (BCArray[cell1]==1)
			       BCArray[cell1]=3;  //mix implicit/explicit boundary
			   }
		       }
		   }
		 else if(_bcMap[fg.id]->bcType == "temperature")
		   {
		     const StorageSite& faces = fg.site;
		     const int faceCount=faces.getCount();
		     const int offSet=faces.getOffset();
		     
		     for(int i=offSet;i<offSet+faceCount;i++)
		       BCfArray[i]=1;
		   }
		 else if(_bcMap[fg.id]->bcType == "Interface")
		   {		 
		     StorageSite& faces0 = fg.site;
		     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces0);
		     const int faceCount=faces0.getCount();
		     const int offSet=faces0.getOffset();
		     
		     for(int i=offSet;i<offSet+faceCount;i++)
		       {
			 BCfArray[i]=4;  //Interface boundary
			 int cell0=BfaceCells(i-offSet,0);
			 int cell1=BfaceCells(i-offSet,1);
			 BCArray[cell0]=2;  //always treated explicitly
			 BCArray[cell1]=4; //interface ghost cells need to be labeled
		       }
		     
		     bool doneAlready=false;
		     foreach(const COMETIC<T>* icPtr, _IClist)
		       {
			 if(icPtr->FgID1==fg.id && icPtr->MeshID1==mesh.getID())
			   {
			     doneAlready=true;
			     break;
			   }
		       }
		     
		     bool foundMatch=false;
		     if(!doneAlready)
		       {
			 StorageSite* faces1Ptr=NULL;
			 int otherFgID,otherMid;
			 for(int otherMeshID=n+1;otherMeshID<numMeshes;otherMeshID++)
			   {
			     const Mesh& otherMesh=*_meshes[otherMeshID];
			     foreach(const FaceGroupPtr otherfgPtr, otherMesh.getBoundaryFaceGroups())
			       {
				 foundMatch=mesh.COMETfindCommonFaces(faces0, otherfgPtr->site, _geomFields);
				 if(foundMatch)
				   {
				     otherFgID=otherfgPtr->id;
				     otherMid=otherMeshID;
				     faces1Ptr=&(otherfgPtr->site);
				     break;
				   }
			       }
			     if(foundMatch)
			       break;
			   }
			 
			 if(foundMatch)
			   {
			     StorageSite& faces1=*faces1Ptr;
			     const Mesh& mesh1=*_meshes[otherMid];
			     COMETIC<T>* icPtr(new COMETIC<T>(n,fg.id,mesh1.getID(),
							      otherFgID, faces1.getCount()));
			     
			     if(_bcMap[fg.id]->InterfaceModel!="DMM")
			       icPtr->InterfaceModel=_bcMap[fg.id]->InterfaceModel;
			     
			     _IClist.push_back(icPtr);
			   }
			 else if(!doneAlready && !foundMatch)
			   {
			     cout<<"Face Group: "<<fg.id<<" MeshID: "<<mesh.getID()<<endl;
			     throw CException("Could not find a matching face group!");
			   }
		       }
		   }// end if interface 
	       }
	   }//end facegroup loop
	 
       }//end meshes loop
     
     //Make map from mesh to IC
     IntArray ICcount(numMeshes);
     ICcount.zero();
     
     foreach(const COMETIC<T>* icPtr, _IClist)
       {
	 ICcount[icPtr->MeshID0]++;
	 ICcount[icPtr->MeshID1]++;
       }

     for(int i=0;i<numMeshes;i++)
       {
	 IntArrPtr MeshToIC(new IntArray(ICcount[i]));
	 _MeshToIC.push_back(MeshToIC);
       }
     
     ICcount.zero();
     const int listSize=_IClist.size();
     for(int ic=0;ic<listSize;ic++)
       {
	 const COMETIC<T>* icPtr=_IClist[ic];
	 (*(_MeshToIC[icPtr->MeshID0]))[ICcount[icPtr->MeshID0]]=ic;
	 (*(_MeshToIC[icPtr->MeshID1]))[ICcount[icPtr->MeshID1]]=ic;
	 ICcount[icPtr->MeshID0]++;
	 ICcount[icPtr->MeshID1]++;
       }

     COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);

     for(int ic=0;ic<listSize;ic++)
       {
	 COMETIC<T>* icPtr=_IClist[ic];
	 if(icPtr->InterfaceModel=="DMM")
	   ComInt.makeDMMcoeffs(*icPtr);
       }
     
     initializeTemperatureBoundaries();
     _residual=updateResid(false);
     cout<<"Creating Coarse Levels..."<<endl;
     MakeCoarseModel(this);
     cout<<"Coarse Levels Completed."<<endl;
  }

  void initCoarse()
  {
    const int numMeshes=_meshes.size();
    const T Tinit=_options["initialTemperature"];

     for (int n=0;n<numMeshes;n++)
       {
	 Mesh& mesh=*_meshes[n];
	 Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	 const int numK=kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 const int kcount=kspace.gettotmodes();

	 shared_ptr<TArray> eArray(new TArray(numcells*kcount));
	 shared_ptr<TArray> e0Array(new TArray(numcells*kcount));
	 shared_ptr<TArray> ResidArray(new TArray(numcells*kcount));
	 shared_ptr<TArray> InjArray(new TArray(numcells*kcount));
	 shared_ptr<TArray> FASArray(new TArray(numcells*kcount));
	 shared_ptr<TArray> tauArray(new TArray(numcells*kcount));
	 kspace.seteArray(eArray);
	 kspace.sete0Array(e0Array);
	 kspace.setResArray(ResidArray);
	 kspace.setInjArray(InjArray);
	 kspace.setFASArray(FASArray);
	 kspace.setTauArray(tauArray);
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 shared_ptr<TArray> deltaTcell(new TArray(numcells));
	 *deltaTcell=0.; 
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 _macro.deltaT.addArray(cells,deltaTcell);

	 shared_ptr<IntArray> f2c(new IntArray(numcells));
	 *f2c=-1;
	 _geomFields.fineToCoarse.addArray(cells, f2c);
	 
	 for(int c=0;c<numcells;c++)
	   {
	     int cellIndex=kspace.getGlobalIndex(c,0);
	     for (int k=0;k<numK;k++)
	       {
		 Tkvol& kv=kspace.getkvol(k);
		 const int numM=kv.getmodenum();

		 for (int m=0;m<numM;m++)
		   {
		     Tmode& mode=kv.getmode(m);
		     const T einit=mode.calce0(Tinit);
		     T tau=mode.gettau();
		     (*eArray)[cellIndex]=einit;
		     (*e0Array)[cellIndex]=einit;
		     (*ResidArray)[cellIndex]=0.;
		     (*InjArray)[cellIndex]=0.;
		     (*FASArray)[cellIndex]=0.;
		     (*tauArray)[cellIndex]=tau;
		     cellIndex++;

		   }
	       }
	     kspace.updateTau(c,Tinit);
	   }

	 BCfaceArray& BCfArray=*(_BFaces[n]);
	 BCcellArray& BCArray=*(_BCells[n]);

	 foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     const StorageSite& faces = fg.site;
	     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);

	     const int faceCount=faces.getCount();
	     const int offSet=faces.getOffset();

	     for(int i=offSet;i<offSet+faceCount;i++)
	       BCfArray[i]=-1;  
	   }

	 foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     if(fg.id>0)
	       {
		 if(_bcMap[fg.id]->bcType == "reflecting")
		   {
		     const StorageSite& faces = fg.site;
		     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces);
		     const int faceCount=faces.getCount();
		     const int offSet=faces.getOffset();
		     const bool Imp=(*(_bcMap[fg.id]))["FullyImplicit"];
		     
		     if(Imp)
		       {
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=2;  //implicit boundary
		       }
		     else
		       {
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=3;  //explicit boundary
		       }
		     
		     if(Imp)
		       {
			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     if(BCArray[cell1]==0)
			       BCArray[cell1]=1;    //implicit boundary only
			     else if(BCArray[cell1]==2)
			       BCArray[cell1]=3;  //mix implicit/explicit boundary
			   }
		       }
		     else
		       {
			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     if(BCArray[cell1]==0)
			       BCArray[cell1]=2;  //explicit boundary only
			     else if (BCArray[cell1]==1)
			       BCArray[cell1]=3;  //mix implicit/explicit boundary
			   }
		       }
		   }
		 else if(_bcMap[fg.id]->bcType == "temperature")
		   {
		     const StorageSite& faces = fg.site;
		     const int faceCount=faces.getCount();
		     const int offSet=faces.getOffset();
		     
		     for(int i=offSet;i<offSet+faceCount;i++)
		       BCfArray[i]=1;
		   }
		 else if(_bcMap[fg.id]->bcType == "Interface")
		   {		 
		     StorageSite& faces0 = fg.site;
		     const CRConnectivity& BfaceCells=mesh.getFaceCells(faces0);
		     const int faceCount=faces0.getCount();
		     const int offSet=faces0.getOffset();
		     
		     for(int i=offSet;i<offSet+faceCount;i++)
		       {
			 BCfArray[i]=4;  //Interface boundary
			 int cell0=BfaceCells(i-offSet,0);
			 BCArray[cell0]=2;  //always treated explicitly
		       }



		     /*
		     bool doneAlready=false;
		     foreach(const COMETIC<T>* icPtr, _IClist)
		       {
			 if(icPtr->FgID1==fg.id && icPtr->MeshID1==mesh.getID())
			   {
			     doneAlready=true;
			     break;
			   }
		       }
		     
		     if(!doneAlready)
		       {
			 bool foundMatch=false;
			 StorageSite* faces1Ptr=NULL;
			 int otherFgID,otherMid;
			 for(int otherMeshID=n+1;otherMeshID<numMeshes;otherMeshID++)
			   {
			     const Mesh& otherMesh=*_meshes[otherMeshID];
			     foreach(const FaceGroupPtr otherfgPtr, otherMesh.getAllFaceGroups())
			       {
				 foundMatch=mesh.COMETfindCommonFaces(faces0, otherfgPtr->site, _geomFields);
				 if(foundMatch)
				   {
				     otherFgID=otherfgPtr->id;
				     otherMid=otherMeshID;
				     faces1Ptr=&(otherfgPtr->site);
				     break;
				   }
			       }
			     if(foundMatch)
			       break;
			   }
			 
			 if(foundMatch)
			   {
			     StorageSite& faces1=*faces1Ptr;
			     const Mesh& mesh1=*_meshes[otherMid];
			     COMETIC<T>* icPtr(new COMETIC<T>(n,fg.id,mesh1.getID(),
							      otherFgID, faces1.getCount()));
			     
			     if(_bcMap[fg.id]->InterfaceModel!="DMM")
			       icPtr->InterfaceModel=_bcMap[fg.id]->InterfaceModel;
			     
			     setLocalScatterMaps(mesh, faces0, mesh1, faces1);
			     
			     _IClist.push_back(icPtr);
			   }
			 else
			   {
			     cout<<"Face Group: "<<fg.id<<" MeshID: "<<mesh.getID()<<endl;
			     throw CException("Could not find a matching face group!");
			   }
		       }
*/
		   }
	       }
	   }
	 
	 TArrptr TlResidCell(new TArray(numcells));
	 *TlResidCell=0.;
	 _macro.TlResidual.addArray(cells,TlResidCell);
	 TArrptr TlInj(new TArray(numcells));
	 *TlInj=0.;
	 _macro.TlInjected.addArray(cells,TlInj);
	 TArrptr TlFAS(new TArray(numcells));
	 *TlFAS=0.;
	 _macro.TlFASCorrection.addArray(cells,TlFAS);
       }

     /*
     const int listSize=_IClist.size();
     for(int ic=0;ic<listSize;ic++)
       {
	 //finer level interface condition
	 COMETIC<T>* icPtr=_IClist[ic];
	 const int meshID0=icPtr->MeshID0;
	 const int meshID1=icPtr->MeshID1;
	 const int fgid0=icPtr->FgID0;
	 const int fgid1=icPtr->FgID1;

	 //this level's interface condition
	 Mesh& mesh0this=*_meshes[meshID0];
	 const FaceGroup& fg0=mesh0this.getFaceGroup(fgid0);
	 const StorageSite& faces0=fg0.site;
	 const int faceCount=faces0.getCount();
	 COMETIC<T>* icPtr2(new COMETIC<T>(meshID0,fgid0,meshID1,
					   fgid1, faceCount));
	 _IClist[ic]=icPtr2;

       }

     
     //Make map from mesh to IC
     IntArray ICcount(numMeshes);
     ICcount.zero();
     
     foreach(const COMETIC<T>* icPtr, _IClist)
       {
	 ICcount[icPtr->MeshID0]++;
	 ICcount[icPtr->MeshID1]++;
       }
     
     for(int i=0;i<numMeshes;i++)
       {
	 IntArrPtr MeshToIC(new IntArray(ICcount[i]));
	 _MeshToIC.push_back(MeshToIC);
       }
     
     ICcount.zero();
     */
     const int listSize=_IClist.size();
     for(int i=0;i<listSize;i++)
       {
	 COMETIC<T> ic=*_IClist[i];
	 const int m0=ic.MeshID0;
	 const int m1=ic.MeshID1;
	 const int fgid0=ic.FgID0;
	 const int fgid1=ic.FgID1;

	 Mesh& mesh0=*_meshes[m0];
	 Mesh& mesh1=*_meshes[m1];

	 FaceGroup& fg0=const_cast<FaceGroup&>(mesh0.getFaceGroup(fgid0));
	 FaceGroup& fg1=const_cast<FaceGroup&>(mesh1.getFaceGroup(fgid1));

	 setLocalScatterMaps(mesh0,fg0.site,mesh1,fg1.site);
       }
     
     /*
     COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);
	 
     for(int ic=0;ic<listSize;ic++)
       {
	 COMETIC<T>* icPtr=_IClist[ic];
	 if(icPtr->InterfaceModel=="DMM")
	   ComInt.makeDMMcoeffs(*icPtr);

       }
     */
     
     applyTemperatureBoundaries();
  }

  void setLocalScatterMaps(const Mesh& mesh0, StorageSite& faces0, const Mesh& mesh1,StorageSite& faces1)
  {//makes the IntArray for the indices of the cells1 to which the cells0 scatter
    
    const IntArray& common01=*(faces0.getCommonMap()[&faces1]);
    IntArrPtr scatter01(new IntArray(faces0.getCount()));
    IntArrPtr scatter10(new IntArray(faces0.getCount()));
    const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
    const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);

    for(int i=0;i<faces0.getCount();i++)
      {
	const int f01=common01[i];
	const int c01=faceCells1(f01,1);
	const int c10=faceCells0(i,1);
	(*scatter01)[i]=c01;
	(*scatter10)[f01]=c10;
      }

    faces0.getScatterMap()[&faces1]=scatter01;
    faces1.getScatterMap()[&faces0]=scatter10;
    
  }

  void applyTemperatureBoundaries()
  {
    const int numMeshes=_meshes.size();

    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*_meshes[n];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if(fg.id>0)
	      {   
		const StorageSite& faces = fg.site;
		const COMETBC<T>& bc = *_bcMap[fg.id];
		
		COMETBoundary<T> cbc(faces, mesh,_geomFields,kspace,_options,fg.id);
		
		if(bc.bcType=="temperature")
		  {	      
		    FloatValEvaluator<T>
		      bTemperature(bc.getVal("specifiedTemperature"),faces);
		    
		    if(_level>0 || _options.Convection=="FirstOrder")
		      cbc.applyTemperatureWallCoarse(bTemperature);
		    else
		      cbc.applyTemperatureWallFine(bTemperature);

		  }
	      }
	  }
      }
  }

  void initializeTemperatureBoundaries()
  {
    const int numMeshes=_meshes.size();

    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*_meshes[n];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if(fg.id>0)
	      {   
		const StorageSite& faces = fg.site;
		const COMETBC<T>& bc = *_bcMap[fg.id];
		
		COMETBoundary<T> cbc(faces, mesh,_geomFields,kspace,_options,fg.id);
		
		if(bc.bcType=="temperature")
		  {	      
		    FloatValEvaluator<T>
		      bTemperature(bc.getVal("specifiedTemperature"),faces);
		    

		    cbc.applyTemperatureWallCoarse(bTemperature);
		    
		    if(_options.Convection=="SecondOrder")
		      cbc.applyTemperatureWallFine(bTemperature);

		  }
	      }
	  }
      }
  }

  void MakeCoarseModel(TCOMET* finerModel)
  {

    if(_options.AgglomerationMethod=="FaceArea")
      {
	int maxLevs=finerModel->getOptions().maxLevels;
	int thisLevel=(finerModel->getLevel())+1;
	
	if(thisLevel<maxLevs)  //assumes # of levels will always work for the mesh
	  {
	    MeshList* newMeshesPtr=new MeshList;
	    TkspList* newKspacesPtr=new TkspList;
	    PhononMacro* newMacroPtr=new PhononMacro("coarse");

	    MakeNewKspaces(finerModel->getKspaces(),*newKspacesPtr);
	    
	    IntArray CoarseCounts(_meshes.size());
	    IntArray CoarseGhost(_meshes.size());
	    map<const StorageSite*, IntArray*> PreFacePairMap;
	    SiteMap siteMap;

	    MakeInteriorCoarseMesh(finerModel->getMeshList(), finerModel->getGeomFields(),
				   *newMeshesPtr, CoarseCounts, PreFacePairMap, CoarseGhost, siteMap);

#ifdef FVM_PARALLEL
	    _geomFields.fineToCoarse.syncLocal();
	    ScatGathMaps coarseScatterMaps;
	    ScatGathMaps coarseGatherMaps;

	    syncGhostCoarsening(finerModel->getMeshList(), finerModel->getGeomFields(),
				*newMeshesPtr, coarseScatterMaps, coarseGatherMaps, CoarseCounts,
				CoarseGhost);

	    makeCoarseScatGath(finerModel->getMeshList(), siteMap, coarseScatterMaps, coarseGatherMaps);
#endif

	    int newCount=FinishCoarseMesh(finerModel->getMeshList(), finerModel->getGeomFields(),
					  *newMeshesPtr, CoarseCounts, PreFacePairMap, CoarseGhost);

	    TCOMET* newModelPtr=new COMETModel(*newMeshesPtr,thisLevel,
					       finerModel->getGeomFields(),
					       *newKspacesPtr,*newMacroPtr);
	    
	    newModelPtr->setFinerLevel(finerModel);
	    finerModel->setCoarserLevel(newModelPtr);
	    newModelPtr->getOptions()=finerModel->getOptions();
	    newModelPtr->getBCs()=finerModel->getBCs();
	    newModelPtr->getMKMap()=finerModel->getMKMap();
	    newModelPtr->getIClist()=finerModel->getIClist();
	    newModelPtr->getMeshICmap()=finerModel->getMeshICmap();
	    COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);
	    ComInt.makeCoarseCoeffs(finerModel->getIClist(),newModelPtr->getIClist(),*newMeshesPtr);
	    newModelPtr->initCoarse();

	    if(_options.DomainStats=="Loud")
	      {
		cout<<"Level: "<<newModelPtr->getLevel()<<endl<<endl;
		newModelPtr->calcDomainStats();
		cout<<endl;
	      }

	    if(newCount>3)
	      newModelPtr->MakeCoarseModel(newModelPtr);
	    else
	      _options.maxLevels=newModelPtr->getLevel();
	  }
      }
    else if(_options.AgglomerationMethod=="AMG")
      throw CException("Have not implemented AMG agglomeration method.");
    else
      throw CException("Unknown agglomeration method.");
  }
  
  void setFinestLevel(TCOMET* finest)
  {
    _finestLevel=finest;
    if(_coarserLevel!=NULL)
      _coarserLevel->setFinestLevel(finest);
  }
  
  void MakeNewKspaces(TkspList& inList, TkspList& outList)
  {
    const int len=inList.size();
    for(int i=0;i<len;i++)
      {
	TkspPtr newKspacePtr=TkspPtr(new Tkspace());
	newKspacePtr->CopyKspace(*inList[i]);
	inList[i]->setCoarseKspace(newKspacePtr);
	outList.push_back(newKspacePtr);
      }
    
    for(int i=0;i<len;i++)
      inList[i]->giveTransmissions();

  }

  void MakeInteriorCoarseMesh(const MeshList& inMeshes, GeomFields& inGeomFields,
			      MeshList& outMeshes, IntArray& CoarseCounts,
			      map<const StorageSite*, IntArray*>& PreFacePairMap, IntArray& CoarseGhost,
			      SiteMap& siteMap)
  {
    
    const int numMeshes=inMeshes.size();
    CoarseCounts=0;

    //coarsen interfaces first
    
    foreach(COMETIC<T>* icPtr, _IClist)
      {
	COMETIC<T>& ic=*icPtr;
	coarsenInterfaceCells(ic, CoarseCounts, inGeomFields, inMeshes);
      }

    //coarsen interiors
    for(int m=0;m<numMeshes;m++)
      {
	const Mesh& mesh=*inMeshes[m];
	const int dim=mesh.getDimension();
	Mesh* newMeshPtr=new Mesh(dim);
	outMeshes.push_back(newMeshPtr);
	newMeshPtr->setID(m);
	siteMap[&(mesh.getCells())]=&(newMeshPtr->getCells());
	int coarseCount=coarsenInterior(m, mesh, inGeomFields, CoarseCounts[m]);
      }

    foreach(COMETIC<T>* icPtr, _IClist)
      {
	COMETIC<T>& ic=*icPtr;
	const int Mid0=ic.MeshID0;
	const int Mid1=ic.MeshID1;
	const int Fid0=ic.FgID0;
	const int Fid1=ic.FgID1;
	const Mesh& mesh0=*inMeshes[Mid0];
	const Mesh& mesh1=*inMeshes[Mid1];
	const StorageSite& faces0=mesh0.getFaceGroup(Fid0).site;
	const StorageSite& faces1=mesh1.getFaceGroup(Fid1).site;
	const StorageSite::CommonMap& ComMap = faces0.getCommonMap();
	const IntArray& common01=*(ComMap.find(&faces1)->second);
	const StorageSite& cells0=mesh0.getCells();
	const StorageSite& cells1=mesh1.getCells();
	const int f1Offset=faces1.getOffset();
	const StorageSite& allFaces0=mesh0.getFaces();
	const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
	const CRConnectivity& cellCells0=mesh0.getCellCells();
	const CRConnectivity& cellCells1=mesh1.getCellCells();
	const CRPtr cellFaces0Ptr=faceCells0.getTranspose();
	const CRConnectivity& cellFaces0=*cellFaces0Ptr;
	const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);
	const CRPtr cellFaces1Ptr=faceCells1.getTranspose();
	const CRConnectivity& cellFaces1=*cellFaces1Ptr;
	const int faceCount=faces0.getCount();
	IntArray& FineToCoarse0=dynamic_cast<IntArray&>(inGeomFields.fineToCoarse[cells0]);
	IntArray& FineToCoarse1=dynamic_cast<IntArray&>(inGeomFields.fineToCoarse[cells1]);
	const int cCount0=CoarseCounts[Mid0];
	const int cCount1=CoarseCounts[Mid1];

	const CRPtr cellsIntGhost0Ptr=cellFaces0.multiply(faceCells0,true);
	const CRConnectivity& cellsIntGhost0=*cellsIntGhost0Ptr;

	Field& FaceAreaField=inGeomFields.area;
	const VectorT3Array& FArea0=
	  dynamic_cast<const VectorT3Array&>(FaceAreaField[allFaces0]);  //global ordering

	IntArray* FaceFine2CoarsePtr0=new IntArray(faceCount);
	*FaceFine2CoarsePtr0=-1;
	IntArray& FaceFine2Coarse0=*FaceFine2CoarsePtr0;
	ic.FineToCoarse0=FaceFine2CoarsePtr0;

	IntArray* FaceFine2CoarsePtr1=new IntArray(faceCount);
	*FaceFine2CoarsePtr1=-1;
	IntArray& FaceFine2Coarse1=*FaceFine2CoarsePtr1;
	ic.FineToCoarse1=FaceFine2CoarsePtr1;

	StorageSite preCoarse0(cCount0);

	CRPtr preCoarseToFineCells0=CRPtr(new CRConnectivity(preCoarse0,cells0));

	preCoarseToFineCells0->initCount();

	for(int c=0;c<cells0.getSelfCount();c++)
	  preCoarseToFineCells0->addCount(FineToCoarse0[c],1);
	
	preCoarseToFineCells0->finishCount();

	for(int c=0;c<cells0.getSelfCount();c++)
	  preCoarseToFineCells0->add(FineToCoarse0[c],c);

	preCoarseToFineCells0->finishAdd();
	
	//initial pairing, no regard to zero summed area
	int coarseFace(0);
	for(int f=0;f<faceCount;f++)
	  {
	    
	    if(FaceFine2Coarse0[f]==-1)
	      {
		const int c00=faceCells0(f,0);   //interior cells
		const int f01=common01[f];
		const int c01=faceCells1(f01,0);
		const int coarseC00=FineToCoarse0[c00];  //coarse indices
		const int coarseC01=FineToCoarse1[c01];

		const int cC00fineCnt=preCoarseToFineCells0->getCount(coarseC00);
		for(int cC00fine=0;cC00fine<cC00fineCnt;cC00fine++)
		  {
		    const int fineCell=(*preCoarseToFineCells0)(coarseC00,cC00fine);
		    if(fineCell!=c00)
		      {
			const int fineGhostCount=cellsIntGhost0.getCount(fineCell);
			for(int fineGhostNum=0;fineGhostNum<fineGhostCount;fineGhostNum++)
			  {
			    const int fineGhost=cellsIntGhost0(fineCell,fineGhostNum);
			    const int f10=cellFaces0(fineGhost,0);
			    const int f11=common01[f10];
			    const int c11=faceCells1(f11,0);
			    
			    if(FineToCoarse1[c11]==coarseC01)
			      {
				if(FaceFine2Coarse0[f10]==-1 && FaceFine2Coarse0[f]==-1)
				  {
				    FaceFine2Coarse0[f10]=coarseFace;
				    FaceFine2Coarse0[f]=coarseFace;
				    FaceFine2Coarse1[f11]=coarseFace;
				    FaceFine2Coarse1[f01]=coarseFace;
				    coarseFace++;
				  }
				else if(FaceFine2Coarse0[f10]==-1 && FaceFine2Coarse0[f]!=-1)
				  {
				    FaceFine2Coarse0[f10]=FaceFine2Coarse0[f];
				    FaceFine2Coarse1[f11]=FaceFine2Coarse0[f01];
				  }
				else if(FaceFine2Coarse0[f10]!=-1 && FaceFine2Coarse0[f]==-1)
				  {
				    FaceFine2Coarse0[f]=FaceFine2Coarse0[f10];
				    FaceFine2Coarse1[f01]=FaceFine2Coarse1[f11];
				  }
			      }
			    
			  }
		      }
		  }

		if(FaceFine2Coarse0[f]==-1)
		  {
		    FaceFine2Coarse0[f]=coarseFace;
		    FaceFine2Coarse1[f01]=coarseFace;
		    coarseFace++;
		  }

	      }
	    //FaceFine2Coarse[f]=coarseFace;
	    //coarseFace++;
	  }


	
	//must check for zero summed area
	VectorT3Array summedArea(coarseFace);
	summedArea.zero();
	for(int f=0;f<faceCount;f++)
	  {
	    const int offset=faces0.getOffset();
	    summedArea[FaceFine2Coarse0[f]]+=FArea0[f+offset];
	  }

	for(int f=0;f<faceCount;f++)
	  {
	    const int f01=common01[f];
	    const int offset=faces0.getOffset();
	    VectorT3& sumVec=summedArea[FaceFine2Coarse0[f]];
	    const VectorT3& partVec=FArea0[f+offset];
	    T sumMag=sqrt(sumVec[0]*sumVec[0]+sumVec[1]*sumVec[1]+
			  sumVec[2]*sumVec[2]);
	    T partMag=sqrt(partVec[0]*partVec[0]+partVec[1]*partVec[1]+
			   partVec[2]*partVec[2]);

	    ///break up face agglomeration if summed area is zero-ish
	    if(sumMag/partMag<1.0)
	      {
		sumVec-=partVec;
		FaceFine2Coarse0[f]=coarseFace;
		FaceFine2Coarse1[f01]=coarseFace;
		coarseFace++;
	      }
	    
	  }

	PreFacePairMap[&faces0]=FaceFine2CoarsePtr0;
	PreFacePairMap[&faces1]=FaceFine2CoarsePtr1;

      }

    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*_meshes[n];
	
	if(!_IClist.empty())
	  {
	    CoarseCounts[n]=correctSingleNeighbor(n, mesh, inGeomFields, 
						  CoarseCounts[n], PreFacePairMap);
	  }

	IntArray& FineToCoarse=dynamic_cast<IntArray&>(
						       inGeomFields.fineToCoarse[mesh.getCells()]);

	CoarseGhost[n]=CoarseCounts[n];
	int outGhost(0);

	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg=*fgPtr;
	    if(fg.id>0)
	      {
		const CRConnectivity& faceCells=mesh.getFaceCells(fg.site);
		const int faceCount=fg.site.getCount();

		if(_bcMap[fg.id]->bcType == "Interface")
		  {
		    IntArray& FaceFine2Coarse=*PreFacePairMap[&fg.site];
		    
		    int coarseFaces(0);
		    for(int i=0;i<faceCount;i++)
		      {
			const int cghost=faceCells(i,1);
			FineToCoarse[cghost]=FaceFine2Coarse[i]+CoarseGhost[n];
			if(FaceFine2Coarse[i]>coarseFaces)
			  coarseFaces=FaceFine2Coarse[i];
		      }
		    CoarseGhost[n]+=coarseFaces+1;
		    outGhost+=coarseFaces+1;
		  }
		else if(_bcMap[fg.id]->bcType == "temperature" ||
			_bcMap[fg.id]->bcType == "reflecting"	)
		  {
		    for(int i=0;i<faceCount;i++)
		      {
			const int cghost=faceCells(i,1);
			FineToCoarse[cghost]=CoarseGhost[n];
			CoarseGhost[n]++;
			outGhost++;
		      }
		  }
	      }
	  }

      }

  }

  void syncGhostCoarsening(const MeshList& inMeshes, GeomFields& inGeomFields,
			   MeshList& outMeshes, ScatGathMaps& coarseScatMaps,
			   ScatGathMaps& coarseGathMaps, IntArray& coarseSizes, 
			   IntArray& CoarseGhost)
  {

    const int numMeshes=inMeshes.size();
    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*inMeshes[n];

	const StorageSite& inCells=mesh.getCells();
	const StorageSite& site=mesh.getCells();
	const int inCellCount=inCells.getSelfCount();
	const int inCellTotal=inCells.getCount();
      
	Field& FineToCoarseField=inGeomFields.fineToCoarse;
	IntArray& coarseIndex=dynamic_cast<IntArray&>(FineToCoarseField[inCells]);
	IntArray tempIndex(inCells.getCountLevel1());
	for(int c=0;c<inCells.getCountLevel1();c++)
	  tempIndex[c]=coarseIndex[c];

	int coarseGhostSize=0;
	int tempGhostSize=0;
	int coarseSize= CoarseGhost[n];

	const StorageSite::GatherMap&  gatherMap  = site.getGatherMap();

	// collect all the toIndices arrays for each storage site from
	// gatherMap
      
	typedef map<const StorageSite*, vector<const Array<int>* > > IndicesMap;
	IndicesMap toIndicesMap;
	//IndicesMap tempIndicesMap;

	foreach(const StorageSite::GatherMap::value_type pos, gatherMap)
	  {
	    const StorageSite& oSite = *pos.first;
	    //const Array<int>& tempIndices = *pos.second;
	    const Array<int>& toIndices = *pos.second;

	    //tempIndicesMap[&oSite].push_back(&tempIndices);
	    toIndicesMap[&oSite].push_back(&toIndices);
	  }
	
	/*
	foreach(IndicesMap::value_type pos, tempIndicesMap)
	  {
	    const StorageSite& oSite = *pos.first;
	    const vector<const Array<int>* > tempIndicesArrays = pos.second;

	    map<int,int> otherToMyMapping;
	    UnorderedSet  gatherSet;

	    foreach(const Array<int>* tempIndicesPtr, tempIndicesArrays)
	      {
		const Array<int>& tempIndices = *tempIndicesPtr;
		const int nGhostRows = tempIndices.getLength();
		for(int ng=0; ng<nGhostRows; ng++)
		  {
		    const int fineIndex = tempIndices[ng];
		    const int coarseOtherIndex = tempIndex[fineIndex];

		    if (coarseOtherIndex < 0)
		      continue;

		    if (otherToMyMapping.find(coarseOtherIndex) !=
			otherToMyMapping.end())
		      {
                      
			tempIndex[fineIndex] = otherToMyMapping[coarseOtherIndex];
		      }
		    else
		      {
			tempIndex[fineIndex] = tempGhostSize+coarseSize;
			otherToMyMapping[coarseOtherIndex] = tempIndex[fineIndex];
			gatherSet.insert( tempIndex[fineIndex] );
			tempGhostSize++;
		      }
		  }
	      }
	  }*/
      
	foreach(IndicesMap::value_type pos, toIndicesMap)
	  {
	    const StorageSite& oSite = *pos.first;
	    const vector<const Array<int>* > toIndicesArrays = pos.second;

          
	    map<int,int> otherToMyMapping;
	    UnorderedSet  gatherSet;

	    foreach(const Array<int>* toIndicesPtr, toIndicesArrays)
	      {
		const Array<int>& toIndices = *toIndicesPtr;
		const int nGhostRows = toIndices.getLength();
		for(int ng=0; ng<nGhostRows; ng++)
		  {
		    const int fineIndex = toIndices[ng];
		    const int coarseOtherIndex = coarseIndex[fineIndex];
                  
		    if (coarseOtherIndex < 0)
		      continue;
                  
		    if (otherToMyMapping.find(coarseOtherIndex) !=
			otherToMyMapping.end())
		      {
			coarseIndex[fineIndex] = otherToMyMapping[coarseOtherIndex];
		      }
		    else
		      {
			coarseIndex[fineIndex] = CoarseGhost[n];
			otherToMyMapping[coarseOtherIndex] = coarseIndex[fineIndex];
			gatherSet.insert( coarseIndex[fineIndex] );
			CoarseGhost[n]++;
		      }
		  
		  }

	      }

	    const int coarseMappersSize = otherToMyMapping.size();

	    shared_ptr<Array<int> > coarseToIndices(new Array<int>(coarseMappersSize));

	    for(int n = 0; n < gatherSet.size(); n++)
	      (*coarseToIndices)[n]   = gatherSet.getData().at(n);
	    
	    SSPair sskey(&site,&oSite);
	    coarseGathMaps [sskey] = coarseToIndices;

	  }
	
	const StorageSite::ScatterMap& scatterMap = site.getScatterMap();
      
	IndicesMap fromIndicesMap;

	foreach(const StorageSite::GatherMap::value_type pos, scatterMap)
	  {
	    const StorageSite& oSite = *pos.first;
	    const Array<int>& fromIndices = *pos.second;

	    fromIndicesMap[&oSite].push_back(&fromIndices);
	  }
      
	foreach(IndicesMap::value_type pos, fromIndicesMap)
	  {
	    const StorageSite& oSite = *pos.first;
	    const vector<const Array<int>* > fromIndicesArrays = pos.second;

	    UnorderedSet  scatterSet;

	    foreach(const Array<int>* fromIndicesPtr, fromIndicesArrays)
	      {
		const Array<int>& fromIndices = *fromIndicesPtr;
		const int nGhostRows = fromIndices.getLength();
		for(int ng=0; ng<nGhostRows; ng++)
		  {
		    const int fineIndex = fromIndices[ng];
		    const int coarseOtherIndex = coarseIndex[fineIndex];
		    if (coarseOtherIndex >= 0)
		      scatterSet.insert( coarseOtherIndex );
		  }

	      }

	    const int coarseMappersSize = scatterSet.size();
          
	    shared_ptr<Array<int> > coarseFromIndices(new Array<int>(coarseMappersSize));
          
	    for(int n = 0; n < scatterSet.size(); n++ ) 
	      (*coarseFromIndices)[n] = scatterSet.getData().at(n);
          
	    SSPair sskey(&site,&oSite);
	    coarseScatMaps[sskey] = coarseFromIndices;
	    
	  }
	
      }
  }

  void makeCoarseScatGath(const MeshList& inMeshes, SiteMap& siteMap,
			   ScatGathMaps& coarseScatMaps, ScatGathMaps& coarseGathMaps)
  {
    const int numMeshes =_meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];
	const StorageSite& fineSite = mesh.getCells();
	StorageSite& coarseSite = *siteMap[&fineSite];
	const StorageSite::ScatterMap& fineScatterMap = fineSite.getScatterMap();
	StorageSite::ScatterMap& coarseScatterMap = coarseSite.getScatterMap();

	foreach(const StorageSite::ScatterMap::value_type& pos, fineScatterMap)
	  {
	    const StorageSite& fineOSite = *pos.first;

	    // the ghost site will not have its corresponding coarse
	    // site created yet so we create it here
	    if (siteMap.find(&fineOSite) == siteMap.end())
	      {
		StorageSite* ghostSite=new StorageSite(-1);
		ghostSite->setGatherProcID (fineOSite.getGatherProcID());
		ghostSite->setScatterProcID(fineOSite.getScatterProcID());
		ghostSite->setTag( fineOSite.getTag() );
		siteMap[&fineOSite]=ghostSite;
	      }
   
	    StorageSite* coarseOSite = siteMap[&fineOSite];
		      
	    SSPair sskey(&fineSite,&fineOSite);
	    coarseScatterMap[coarseOSite] = coarseScatMaps[sskey];

	  }
		  
	const StorageSite::GatherMap& fineGatherMap = fineSite.getGatherMap();
	StorageSite::GatherMap& coarseGatherMap = coarseSite.getGatherMap();

	foreach(const StorageSite::GatherMap::value_type& pos, fineGatherMap)
	  {
	    const StorageSite& fineOSite = *pos.first;
	    StorageSite& coarseOSite = *siteMap[&fineOSite];
	    SSPair sskey(&fineSite,&fineOSite);

	    coarseGatherMap[&coarseOSite] = coarseGathMaps[sskey];
	  }
		  
      }
  }

  int FinishCoarseMesh(const MeshList& inMeshes, GeomFields& inGeomFields,
		       MeshList& outMeshes, IntArray& CoarseCounts,
		       map<const StorageSite*, IntArray*>& PreFacePairMap,
		       IntArray& coarseGhost)
  {

    const int numMeshes=inMeshes.size();

    int smallestMesh=-1;
    
    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*inMeshes[n];
	Mesh* newMeshPtr=outMeshes[n];
	int coarseCount=CoarseCounts[n];
	const StorageSite& inCells=mesh.getCells();
	IntArray& FineToCoarse=dynamic_cast<IntArray&>(inGeomFields.fineToCoarse[inCells]);
	StorageSite& outCells=newMeshPtr->getCells();
	StorageSite& outFaces=newMeshPtr->getFaces();
	const StorageSite& inFaces=mesh.getFaces();
	const int inCellTotal=inCells.getCount();
	const int inFaceCount=inFaces.getCount();

	const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
	Field& areaMagField=inGeomFields.areaMag;
	const TArray& areaMagArray=
	  dynamic_cast<const TArray&>(areaMagField[inFaces]);

	Field& FaceAreaField=inGeomFields.area;
	const VectorT3Array& inFA=
	  dynamic_cast<const VectorT3Array&>(FaceAreaField[inFaces]);
	
	//make the coarse cell to fine cell connectivity.
	outCells.setCount(CoarseCounts[n],coarseGhost[n]-CoarseCounts[n]);
	CRPtr CoarseToFineCells=CRPtr(new CRConnectivity(outCells,inCells));
	CoarseToFineCells->initCount();
	
	for(int c=0;c<inCellTotal;c++)
	  CoarseToFineCells->addCount(FineToCoarse[c],1);
	
	CoarseToFineCells->finishCount();

	for(int c=0;c<inCellTotal;c++)
	  CoarseToFineCells->add(FineToCoarse[c],c);

	CoarseToFineCells->finishAdd();

	//connectivity between itself (cells) and its finer mesh cells. 
	newMeshPtr->setConnectivity(outCells,inCells,CoarseToFineCells);

	CRPtr FineFacesCoarseCells=CRPtr(new CRConnectivity(inFaces,outCells));
	FineFacesCoarseCells->initCount();

	//count surviving faces
	for(int f=0;f<inFaceCount;f++)
	  {
	    int coarse0=FineToCoarse[inFaceinCells(f,0)];
	    int coarse1=FineToCoarse[inFaceinCells(f,1)];
	    if(coarse0!=coarse1)
	      FineFacesCoarseCells->addCount(f,2);
	  }
	
	FineFacesCoarseCells->finishCount();
	
	//make non-zero's
	for(int f=0;f<inFaceCount;f++)
	  {
	    int fc0=inFaceinCells(f,0);
	    int fc1=inFaceinCells(f,1);
	    int cc0=FineToCoarse[fc0];
	    int cc1=FineToCoarse[fc1];
	    if(cc0!=cc1)
	      {
		FineFacesCoarseCells->add(f,cc0);
		FineFacesCoarseCells->add(f,cc1);
	      }
	  }

	FineFacesCoarseCells->finishAdd();

	CRPtr CoarseCellsFineFaces=FineFacesCoarseCells->getTranspose();
	CRPtr CellCellCoarse=CoarseCellsFineFaces->multiply(*FineFacesCoarseCells,true);

	int counter=0;  //coarse face counter
	BArray counted(outCells.getCount());
	counted=false;
	for(int c=0;c<outCells.getCount();c++)
	  {
	    counted[c]=true;
	    const int neibs=CellCellCoarse->getCount(c);

	    for(int n=0;n<neibs;n++)
	      {
		const int c1=(*CellCellCoarse)(c,n);
		if(neibs>1)
		  {
		    if(!counted[c1])
		      {
			counter++;
		      }
		  }
		else  //gives two faces to cells with 1 neighbor
		  {
		    if(c>=outCells.getSelfCount())
		      {
			if(!counted[c1])
			  counter++;
		      }
		    else  //interior cells
		      {
			if(counted[c1])
			  counter++;
			else
			  counter+=2;
		      }
		  }
	      }
	  }

	outFaces.setCount(counter);

	CRPtr CoarseCellCoarseFace=CRPtr(new CRConnectivity(outCells,outFaces));
	CoarseCellCoarseFace->initCount();

	int tempCount(0);
	for(int c=0;c<outCells.getCount();c++)
	  {
	    const int neibs=CellCellCoarse->getCount(c);
	    CoarseCellCoarseFace->addCount(c,neibs);
	    tempCount+=neibs;
	    if((neibs==1) && (c<outCells.getSelfCount()))   //adding 2 faces to cells with one neighbor
	      {
		CoarseCellCoarseFace->addCount(c,neibs);
		tempCount+=neibs;
	      }
	  }

	CoarseCellCoarseFace->finishCount();

	//make cell connectivity to interior faces.
	counter=0;
	counted=false;
	for(int c=0;c<outCells.getSelfCount();c++)
	  {
	    counted[c]=true;
	    const int neibs=CellCellCoarse->getCount(c);
	    for(int n=0;n<neibs;n++)
	      {
		const int c1=(*CellCellCoarse)(c,n);
		if(neibs>1)
		  {
		    if(!counted[c1] && c1<outCells.getSelfCount())
		      {
			int increment=CoarseCellCoarseFace->add(c,counter);
			increment=CoarseCellCoarseFace->add(c1,counter);
			counter++;
		      }
		  }
		else
		  {
		    if(counted[c1] && c1<outCells.getSelfCount())
		      {
			CoarseCellCoarseFace->add(c,counter);
			CoarseCellCoarseFace->add(c1,counter);
			counter++;
		      }
		    else if(!counted[c1] && c1<outCells.getSelfCount())
		      {
			CoarseCellCoarseFace->add(c,counter);
			CoarseCellCoarseFace->add(c1,counter);
			counter++;
			CoarseCellCoarseFace->add(c,counter);
			CoarseCellCoarseFace->add(c1,counter);
			counter++;
		      }
		  }
	      }
	  }

	const int coarseInteriorCount=counter;

	//cell connectivity to boundary faces
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg=*fgPtr;
	    const int faceCount=fg.site.getCount();
	    const CRConnectivity& inBfaceCells=mesh.getFaceCells(fg.site);

	    if(fg.id>0)
	      {
		int faceCount(0);
		if(_bcMap[fg.id]->bcType == "Interface")
		  {
		    IntArray& FaceFine2Coarse=*PreFacePairMap[&fg.site];
		    int coarseFaces(0);
		    for(int i=0;i<fg.site.getCount();i++)
		      if(FaceFine2Coarse[i]>coarseFaces)
			coarseFaces=FaceFine2Coarse[i];
		    faceCount=coarseFaces+1;

		    BArray countedFaces(faceCount);
		    countedFaces=false;
		    for(int i=0;i<fg.site.getCount();i++)
		      {
			const int coarseFace=FaceFine2Coarse[i];
			if(!countedFaces[coarseFace])
			  {
			    const int globFace=i+fg.site.getOffset();
			    const int cc1=(*FineFacesCoarseCells)(globFace,0);
			    const int cc2=(*FineFacesCoarseCells)(globFace,1);
			    int increment=CoarseCellCoarseFace->add(cc1,counter);
			    increment=CoarseCellCoarseFace->add(cc2,counter);
			    counter++;
			    countedFaces[coarseFace]=true;
			  }
		      }

		  }
		else
		  {
		    faceCount=fg.site.getCount();
		    for(int i=0;i<faceCount;i++)
		      {
			const int c1=inBfaceCells(i,1);   //fine ghost cell
			const int cc1=FineToCoarse[c1];   //coarse ghost cell
			const int cc2=(*CellCellCoarse)(cc1,0);  //coarse interior cell
			int increment=CoarseCellCoarseFace->add(cc1,counter);
			increment=CoarseCellCoarseFace->add(cc2,counter);
			counter++;
		      }
		  }
	      }
	  }
	
	// have to count coarse interfaces (parallel) -- !!

	foreach(const StorageSite::GatherMap::value_type pos, outCells.getGatherMap())
	  {
	    const StorageSite& oSite = *pos.first;
	    const Array<int>& toIndices = *pos.second;

	    int to_where  = oSite.getGatherProcID();
	    if ( to_where != -1 )
	      {
		const int fgCount=toIndices.getLength();
		
		for(int i=0;i<fgCount;i++)
		  {
		    const int c1=toIndices[i];
		    const int c1neibs=CellCellCoarse->getCount(c1);
		    for(int j=0;j<c1neibs;j++)
		      {
			const int c2=(*CellCellCoarse)(toIndices[i],j);
			CoarseCellCoarseFace->add(c2,counter);
			CoarseCellCoarseFace->add(toIndices[i],counter);
			counter++;
		      }
		  }
	      }
	  }

	CoarseCellCoarseFace->finishAdd();

	CRPtr CoarseFaceCoarseCell=CoarseCellCoarseFace->getTranspose();

	newMeshPtr->setConnectivity(outCells,outFaces,CoarseCellCoarseFace);
	newMeshPtr->setConnectivity(outFaces,outCells,CoarseFaceCoarseCell);

	CRPtr CoarseFacesFineFaces=CRPtr(new CRConnectivity(outFaces,inFaces));
	CoarseFacesFineFaces->initCount();

	VectorT3Array areaSum(outFaces.getCount());
	areaSum.zero();

	for(int f=0;f<inFaceCount;f++)
	  {
	    int fc0=inFaceinCells(f,0);
	    int fc1=inFaceinCells(f,1);
	    const int cc0=FineToCoarse[fc0];
	    const int cc1=FineToCoarse[fc1];

	    int cc0neibs=CellCellCoarse->getCount(cc0);
	    int cc1neibs=CellCellCoarse->getCount(cc1);

	    if(cc0>=outCells.getSelfCount())
	      cc0neibs++;
	    if(cc1>=outCells.getSelfCount())
	      cc1neibs++;

	    if(cc1!=cc0)
	      {
		const int cfaces=CoarseCellCoarseFace->getCount(cc0);
		
		for(int cf=0;cf<cfaces;cf++)
		  {
		    const int face=(*CoarseCellCoarseFace)(cc0,cf);
		    const int tempc0=(*CoarseFaceCoarseCell)(face,0);
		    const int tempc1=(*CoarseFaceCoarseCell)(face,1);
		    
		    if(((cc0==tempc0)&&(cc1==tempc1))||((cc1==tempc0)&&(cc0==tempc1)))
		      {
			T sign(0);

			if(cc1==tempc0)
			   sign*=-1.;

			areaSum[face]+=inFA[f]*sign;

			if(cc1neibs>1 && cc0neibs>1)
			  {
			    CoarseFacesFineFaces->addCount(face,1);
			    break;
			  }
			else
			  {
			    T sumMag=sqrt(areaSum[face][0]*areaSum[face][0]+
					  areaSum[face][1]*areaSum[face][1]+
					  areaSum[face][2]*areaSum[face][2]);

			    T partMag=sqrt(inFA[f][0]*inFA[f][0]+
					   inFA[f][1]*inFA[f][1]+
					   inFA[f][2]*inFA[f][2]);
			
			    if(sumMag/partMag>.75)
			      {
				CoarseFacesFineFaces->addCount(face,1);
				break;
			      }
			    else
			      areaSum[face]-=inFA[f]*sign;

			  }
		      }
		  }
	      }
	  }

	CoarseFacesFineFaces->finishCount();

	areaSum.zero();

	for(int f=0;f<inFaceCount;f++)
	  {
	    int fc0=inFaceinCells(f,0);
	    int fc1=inFaceinCells(f,1);
	    const int cc0=FineToCoarse[fc0];
	    const int cc1=FineToCoarse[fc1];

	    int cc0neibs=CellCellCoarse->getCount(cc0);
	    int cc1neibs=CellCellCoarse->getCount(cc1);

	    if(cc0>=outCells.getSelfCount())
	      cc0neibs++;
	    if(cc1>=outCells.getSelfCount())
	      cc1neibs++;
	    
	    if(cc1!=cc0)
	      {
		const int cfaces=CoarseCellCoarseFace->getCount(cc0);
		
		for(int cf=0;cf<cfaces;cf++)
		  {
		    const int face=(*CoarseCellCoarseFace)(cc0,cf);
		    const int tempc0=(*CoarseFaceCoarseCell)(face,0);
		    const int tempc1=(*CoarseFaceCoarseCell)(face,1);
		    
		    if(((cc0==tempc0)&&(cc1==tempc1))||((cc1==tempc0)&&(cc0==tempc1)))
		      {

			T sign(1);
			if(cc1==tempc0)
			  sign*=-1.;

			areaSum[face]+=inFA[f]*sign;

			if(cc1neibs>1 && cc0neibs>1)
			  {
			    CoarseFacesFineFaces->add(face,f);
			    break;
			  }
			else
			  {
			    T sumMag=sqrt(areaSum[face][0]*areaSum[face][0]+
					  areaSum[face][1]*areaSum[face][1]+
					  areaSum[face][2]*areaSum[face][2]);
			    
			    T partMag=sqrt(inFA[f][0]*inFA[f][0]+
					   inFA[f][1]*inFA[f][1]+
					   inFA[f][2]*inFA[f][2]);
			    
			    if(sumMag/partMag>.75)
			      {
				CoarseFacesFineFaces->add(face,f);
				break;
			      }
			    else
			      areaSum[face]-=inFA[f]*sign;
			  }
		      }
		  }
	      }
	  }

	CoarseFacesFineFaces->finishAdd();
	
	//const StorageSite& interiorFaces=
	newMeshPtr->createInteriorFaceGroup(coarseInteriorCount);

	int inOffset=coarseInteriorCount;
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg=*fgPtr;
	    const int faceCount=fg.site.getCount();
	    if(fg.id>0)
	      {
		if(_bcMap[fg.id]->bcType == "Interface")
		  {
		    IntArray& FaceFine2Coarse=*PreFacePairMap[&fg.site];
		    int coarseFaces(0);
		    for(int i=0;i<faceCount;i++)
		      if(FaceFine2Coarse[i]>coarseFaces)
			coarseFaces=FaceFine2Coarse[i];

		    coarseFaces+=1;

		    //const StorageSite& newBoundarySite=
		    newMeshPtr->createBoundaryFaceGroup(coarseFaces,inOffset,fg.id,fg.groupType);
		    inOffset+=coarseFaces;
		  }
		else
		  {
		    const int size=fg.site.getCount();
		    const StorageSite& newBoundarySite=newMeshPtr->createBoundaryFaceGroup(
											   size,inOffset,fg.id,fg.groupType);
		    const VectorT3Array& oldFgCoords=dynamic_cast<const VectorT3Array&>(inGeomFields.coordinate[fg.site]);
		    VT3Ptr newFgCoordsPtr=VT3Ptr(new VectorT3Array(size));
		    (*newFgCoordsPtr)=oldFgCoords;
		    inGeomFields.coordinate.addArray(newBoundarySite,newFgCoordsPtr);
		    inOffset+=size;
		  }
	      }
	  }

	// have to make coarse interface face groups -- !!

	const StorageSite::GatherMap&  coarseGatherMap  = outCells.getGatherMap();

	foreach(const StorageSite::GatherMap::value_type pos, coarseGatherMap)
	  {
	    const StorageSite& oSite = *pos.first;
	    const Array<int>& toIndices = *pos.second;

	    int to_where  = oSite.getGatherProcID();
	    if ( to_where != -1 )
	      {
		const int cellCount=toIndices.getLength();
		int faceCount(0);

		for(int i=0;i<cellCount;i++)
		  {
		    const int c1=toIndices[i];
		    faceCount+=CellCellCoarse->getCount(c1);
		  }

		const StorageSite& newBoundarySite=
		  newMeshPtr->createInterfaceGroup(faceCount,inOffset,to_where);
		inOffset+=faceCount;
	      }
	  }
	
	//now make the geom fields
	const int outCellsCount=outCells.getSelfCount();
	const int outCellsTotCount=outCells.getCount();
	VT3Ptr outCellCoordsPtr=VT3Ptr(new VectorT3Array(outCellsTotCount));
	TArrptr outCellVolumePtr=TArrptr(new TArray(outCellsTotCount));
	TArray& outCV=*outCellVolumePtr;
	VectorT3Array& outCoords=*outCellCoordsPtr;
	outCV=0.;

	Field& VolumeField=inGeomFields.volume;
	const TArray& inCV=dynamic_cast<const TArray&>(VolumeField[inCells]);
	const VectorT3Array& inCoords=dynamic_cast<const VectorT3Array&>(inGeomFields.coordinate[inCells]);

	for(int c=0;c<outCellsCount;c++)
	  {
	    const int fineCount=CoarseToFineCells->getCount(c);
	    VectorT3 newCoord;
	    newCoord[0]=0.;
	    newCoord[1]=0.;
	    newCoord[2]=0.;
	    for(int i=0;i<fineCount;i++)
	      {
		int fc=(*CoarseToFineCells)(c,i);
		outCV[c]+=inCV[fc];
		newCoord+=inCoords[fc]*inCV[fc];
	      }
	    outCoords[c]=newCoord/outCV[c];
	  }

	for(int c=outCellsCount;c<outCellsTotCount;c++)
	  {
	    const int fineCount=CoarseToFineCells->getCount(c);
	    VectorT3 newCoord;
	    newCoord[0]=0.;
	    newCoord[1]=0.;
	    newCoord[2]=0.;
	    for(int i=0;i<fineCount;i++)
	      {
		int fc=(*CoarseToFineCells)(c,i);
		newCoord+=inCoords[fc];
	      }
	    outCoords[c]=newCoord;  ///needs to be changed
	  }	
	
	VolumeField.addArray(outCells,outCellVolumePtr);
	inGeomFields.coordinate.addArray(outCells,outCellCoordsPtr);

	const int outFacesCount=outFaces.getCount();
	VT3Ptr outFaceAreaPtr=VT3Ptr(new VectorT3Array(outFacesCount));
	VectorT3Array& outFA=*outFaceAreaPtr;
	TArrptr outFaceAreaMagPtr=TArrptr(new TArray(outFacesCount));
	TArray& outFAMag=*outFaceAreaMagPtr;

	VectorT3 myZero;
	myZero[0]=0.;
	myZero[1]=0.;
	myZero[2]=0.;

	outFA=myZero;
	outFAMag=0.;
	for(int f=0;f<outFacesCount;f++)
	  {
	    const int fineCount=CoarseFacesFineFaces->getCount(f);
	    const int cCell0=(*CoarseFaceCoarseCell)(f,0);
	    for(int i=0;i<fineCount;i++)
	      {
		const int fFace=(*CoarseFacesFineFaces)(f,i);
		const int fCell0=inFaceinCells(fFace,0);
		const int CCell0=FineToCoarse[fCell0];
		
		//must make sure the area vector is pointing
		//from c0 to c1
		if(CCell0==cCell0)
		  outFA[f]+=inFA[fFace];
		else
		  outFA[f]-=inFA[fFace];
		outFAMag[f]+=areaMagArray[fFace];
	      }
	  }

	FaceAreaField.addArray(outFaces,outFaceAreaPtr);
	areaMagField.addArray(outFaces,outFaceAreaMagPtr);
	//cout<<"Level: "<<_level+1<<" Mesh: "<<n<<" Cells: "<<outCells.getSelfCount()<<endl;

	if(smallestMesh<0)
	  smallestMesh=outCells.getSelfCount();
	else
	  {
	    if(outCells.getSelfCount()<smallestMesh)
	      smallestMesh=outCells.getSelfCount();
	  }

	/*
	//This is for checking purposes only
	cout<<"Coarse Faces to Fine Faces"<<endl;
	for(int f=0;f<outFaces.getCount();f++)
	  {
	    const int neibs=CoarseFacesFineFaces->getCount(f);
	    for(int n=0;n<neibs;n++)
	      cout<<f<<" "<<(*CoarseFacesFineFaces)(f,n)<<endl;
	    cout<<endl;
	  }
	cout<<"Coarse Cells to Coarse Faces"<<endl;
	for(int c=0;c<outCells.getCount();c++)
	  {
	    const int neibs=CoarseCellCoarseFace->getCount(c);
	    for(int n=0;n<neibs;n++)
	      cout<<c<<" "<<(*CoarseCellCoarseFace)(c,n)<<endl;
	    cout<<endl;
	  }
	*/
      }

    return smallestMesh;
  }

  int coarsenInterior(const int m, const Mesh& mesh, 
		      GeomFields& inGeomFields, int& coarseCount)
  {
    const StorageSite& inCells=mesh.getCells();
    const StorageSite& inFaces=mesh.getFaces();
    const int inCellCount=inCells.getSelfCount();
    const CRConnectivity& inCellinFaces=mesh.getCellFaces();
    const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
    Field& areaMagField=inGeomFields.areaMag;
    const TArray& areaMagArray=dynamic_cast<const TArray&>(areaMagField[inFaces]);
    const BCfaceArray& inBCfArray=*(_BFaces[m]);

    IntArray& FineToCoarse=dynamic_cast<IntArray&>(inGeomFields.fineToCoarse[inCells]);

    //first sweep to make initial pairing
    int pairWith;
    for(int c=0;c<inCellCount;c++)
      {
	if(FineToCoarse[c]<0) //dont bother if im already paired
	  {
	    //loop through all neighbors to find pairing
	    const int neibCount=inCellinFaces.getCount(c);
	    pairWith=-1;
	    T maxArea=0.;
	    int c2;
	    for(int neib=0;neib<neibCount;neib++)
	      {
		const int f=inCellinFaces(c,neib);
		    
		if(inBCfArray[f]==0)  //not a boundary face
		  {
		    if(c==inFaceinCells(f,1))
		      c2=inFaceinCells(f,0);
		    else
		      c2=inFaceinCells(f,1);
			
		    if(FineToCoarse[c2]==-1)  //not already paired
		      {
			if(areaMagArray[f]>maxArea)
			  {
			    pairWith=c2;
			    maxArea=areaMagArray[f];
			  }
		      }

		  }
	      }
		
	    if(pairWith!=-1)
	      {
		FineToCoarse[c]=coarseCount;
		FineToCoarse[pairWith]=coarseCount;
		coarseCount++;
	      }
	  }
      }
	
    //second sweep to group stragglers, or group with self
    for(int c=0;c<inCellCount;c++)
      {
	if(FineToCoarse[c]==-1)
	  {
	    const int neibCount=inCellinFaces.getCount(c);
	    T maxArea=0.;
	    int c2,c2perm;
	    pairWith=-2;
		
	    for(int neib=0;neib<neibCount;neib++)
	      {
		const int f=inCellinFaces(c,neib);
		    
		if(inBCfArray[f]==0)  //not a boundary face
		  {
		    if(c==inFaceinCells(f,1))
		      c2=inFaceinCells(f,0);
		    else
		      c2=inFaceinCells(f,1);
			
		    if(areaMagArray[f]>maxArea)
		      {
			pairWith=FineToCoarse[c2]; //coarse level cell
			c2perm=c2;                 //fine level cell
			maxArea=areaMagArray[f];
		      }
		  }
	      }
		
	    if(pairWith==-2)  //could not find suitable pair
	      {
		FineToCoarse[c]=coarseCount;
		coarseCount++;
	      }
	    else  //found cell to pair with
	      {
		if(FineToCoarse[c2perm]==-1)
		  {
		    FineToCoarse[c]=coarseCount;
		    FineToCoarse[c2perm]=coarseCount;
		    coarseCount++;
		  }
		else
		  FineToCoarse[c]=FineToCoarse[c2perm];
	      }
	  }	    
      }
    return coarseCount;
  }

  int correctSingleNeighbor(const int m, const Mesh& mesh,
			    GeomFields& inGeomFields, int coarseCount,
			    map<const StorageSite*, IntArray*> PreFacePairMap)
  {

    const StorageSite& inCells=mesh.getCells();
    const int inCellTotal=inCells.getCount();
    const StorageSite& inFaces=mesh.getFaces();
    const int inFaceCount=inFaces.getCount();
    const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
    const IntArray& FineToCoarseConst=dynamic_cast<const IntArray&>
      (inGeomFields.fineToCoarse[inCells]);

    IntArray FineToCoarse(FineToCoarseConst.getLength());
    FineToCoarse=FineToCoarseConst;
    
    int coarseGhost=coarseCount;
    int outGhost(0);
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg=*fgPtr;

	if(fg.id>0)
	  {
	    const CRConnectivity& faceCells=mesh.getFaceCells(fg.site);
	    const int faceCount=fg.site.getCount();

	    if(_bcMap[fg.id]->bcType == "Interface")
	      {
		IntArray& FaceFine2Coarse=*PreFacePairMap[&fg.site];
		    
		int coarseFaces(0);
		for(int i=0;i<faceCount;i++)
		  {
		    const int cghost=faceCells(i,1);
		    FineToCoarse[cghost]=FaceFine2Coarse[i]+coarseGhost;
		    if(FaceFine2Coarse[i]>coarseFaces)
		      coarseFaces=FaceFine2Coarse[i];
		  }
		coarseGhost+=coarseFaces+1;
		outGhost+=coarseFaces+1;
	      }
	    else if(_bcMap[fg.id]->bcType == "temperature" ||
		    _bcMap[fg.id]->bcType == "reflecting"	)
	      {
		for(int i=0;i<faceCount;i++)
		  {
		    const int cghost=faceCells(i,1);
		    FineToCoarse[cghost]=coarseGhost;
		    coarseGhost++;
		    outGhost++;
		  }
	      }
	  
	  }
      }

    foreach(const FaceGroupPtr fgPtr, mesh.getInterfaceGroups())
      {
	const FaceGroup& fg=*fgPtr;
	const CRConnectivity& faceCells=mesh.getFaceCells(fg.site);
	const int faceCount=fg.site.getCount();
	
	for(int i=0;i<faceCount;i++)
	  {
	    const int cghost=faceCells(i,1);
	    FineToCoarse[cghost]=coarseGhost;
	    coarseGhost++;
	    outGhost++;
	  }
	
      }

    StorageSite preOutCells(coarseCount,outGhost);
    CRPtr CoarseToFineCells=CRPtr(new CRConnectivity(preOutCells,inCells));
    CoarseToFineCells->initCount();
	
    for(int c=0;c<inCellTotal;c++)
      CoarseToFineCells->addCount(FineToCoarse[c],1);
	
    CoarseToFineCells->finishCount();

    for(int c=0;c<inCellTotal;c++)
      CoarseToFineCells->add(FineToCoarse[c],c);

    CoarseToFineCells->finishAdd();

    CRPtr FineFacesCoarseCells=CRPtr(new CRConnectivity(inFaces,preOutCells));
    FineFacesCoarseCells->initCount();

    //count surviving faces
    int survivingFaces=0;
    int coarse0, coarse1;
    for(int f=0;f<inFaceCount;f++)
      {
	coarse0=FineToCoarse[inFaceinCells(f,0)];
	coarse1=FineToCoarse[inFaceinCells(f,1)];
	if(coarse0!=coarse1)
	  {
	    survivingFaces++;
	    FineFacesCoarseCells->addCount(f,2);
	  }
      }
	
    FineFacesCoarseCells->finishCount();
	
    //make non-zero's
    for(int f=0;f<inFaceCount;f++)
      {
	int fc0=inFaceinCells(f,0);
	int fc1=inFaceinCells(f,1);
	int cc0=FineToCoarse[fc0];
	int cc1=FineToCoarse[fc1];
	if(cc0!=cc1)
	  {
	    FineFacesCoarseCells->add(f,cc0);
	    FineFacesCoarseCells->add(f,cc1);
	  }
      }

    FineFacesCoarseCells->finishAdd();

    CRPtr CoarseCellsFineFaces=FineFacesCoarseCells->getTranspose();
    CRPtr CellCellCoarse=CoarseCellsFineFaces->multiply(*FineFacesCoarseCells,true);


    for(int c=0;c<preOutCells.getSelfCount();c++)
      {
	if(CellCellCoarse->getCount(c)==1)
	  {
	    /*
	    const int bigCoarse=(*CellCellCoarse)(c,0);
	    const int fineCount=CoarseToFineCells->getCount(c);

	    for(int j=0;j<fineCount;j++)
	      {
		const int fc=(*CoarseToFineCells)(c,j);
		FineToCoarse[fc]=bigCoarse;
	      }
	    coarseCount-=1;
	    */
	    
	    const int removal=(*CoarseToFineCells)(c,0);
	    FineToCoarse[removal]=coarseCount;
	    coarseCount++;
	    
	  }
      }

    return coarseCount;
  }

  void coarsenInterfaceCells(COMETIC<T>& ic, IntArray& coarseCounts, 
			     GeomFields& inGeomFields, const MeshList& inMeshes)
  {
    const int Mid0=ic.MeshID0;
    const int Mid1=ic.MeshID1;
    const int Fid0=ic.FgID0;
    const int Fid1=ic.FgID1;
    const Mesh& mesh0=*inMeshes[Mid0];
    const Mesh& mesh1=*inMeshes[Mid1];
    const StorageSite& faces0=mesh0.getFaceGroup(Fid0).site;
    const StorageSite& faces1=mesh1.getFaceGroup(Fid1).site;
    const StorageSite& cells0=mesh0.getCells();
    const StorageSite& cells1=mesh1.getCells();
    const StorageSite& globFaces0=mesh0.getFaces();
    const StorageSite::CommonMap& ComMap = faces0.getCommonMap();
    const IntArray& common01=*(ComMap.find(&faces1)->second);
    const CRConnectivity& faceCells0=mesh0.getFaceCells(faces0);
    const CRConnectivity& faceCells1=mesh1.getFaceCells(faces1);
    const CRConnectivity& cellCells0=mesh0.getCellCells();
    const CRConnectivity& cellCells1=mesh1.getCellCells();
    const CRConnectivity& globFaceCells0=mesh0.getAllFaceCells();
    const CRConnectivity& globCellFaces0=mesh0.getCellFaces();
    const int inFaceCount=faces0.getCount();
    Field& areaMagField=inGeomFields.areaMag;
    const TArray& areaMagArray0=dynamic_cast<const TArray&>(areaMagField[globFaces0]);
    const BCfaceArray& inBCfArray0=*(_BFaces[Mid0]);

    const CRPtr cellFaces0Ptr=faceCells0.getTranspose();
    const CRConnectivity& cellFaces0=*cellFaces0Ptr;
    const CRPtr cellsIntGhost0Ptr=cellFaces0.multiply(faceCells0,true);
    const CRConnectivity& cellsIntGhost0=*cellsIntGhost0Ptr;

    IntArray CoarseFineCount(cells0.getCount());
    CoarseFineCount.zero();

    IntArray& f2c0=dynamic_cast<IntArray&>(inGeomFields.fineToCoarse[cells0]);
    IntArray& f2c1=dynamic_cast<IntArray&>(inGeomFields.fineToCoarse[cells1]);
    int& count0=coarseCounts[Mid0];
    int& count1=coarseCounts[Mid1];

    const VectorT3Array& cpos0=dynamic_cast<const VectorT3Array&>(inGeomFields.coordinate[cells0]);
    const VectorT3Array& cpos1=dynamic_cast<const VectorT3Array&>(inGeomFields.coordinate[cells1]);
    
    //agglomerate mesh0 first
    int pairWith;

    if(_level>=0)
      {
	for(int f=0;f<inFaceCount;f++)
	  {
	    const int c=faceCells0(f,0);
	    const int cCoarse=CoarseFineCount[c];
	    if(f2c0[c]<0 || cCoarse<2) //dont bother if im already paired
	      {
		//loop through all neighbors to find pairing
		const int neibCount=cellCells0.getCount(c);
		pairWith=-2;
		T maxArea=0.;
		int c2;
		for(int neib=0;neib<neibCount;neib++)
		  {
		    const int fglob=globCellFaces0(c,neib);
		    
		    if(inBCfArray0[fglob]==0)  //not a boundary face
		      {
			if(c==globFaceCells0(fglob,1))
			  c2=globFaceCells0(fglob,0);
			else
			  c2=globFaceCells0(fglob,1);
			
			int c2Coarse=f2c0[c2];
			
			if(_level==0)
			  {
			    if(c2Coarse==-1 || CoarseFineCount[c2Coarse]<5)  //not already paired
			      {
				if(areaMagArray0[fglob]>maxArea)
				  {
				    pairWith=c2;
				    maxArea=areaMagArray0[fglob];
				    if(cellsIntGhost0.getCount(c2)>0)  //pick first interface cell
				      break;
				  }
			      }
			  }
			else
			  {
			    if(c2Coarse==-1)
			      {
				if(areaMagArray0[fglob]>maxArea)
				  {
				    pairWith=c2;
				    maxArea=areaMagArray0[fglob];
				    if(cellsIntGhost0.getCount(c2)>0)  //pick first interface cell
				      break;
				  }
			      }
			  }

		      }
		  }
		
		if(pairWith>=0)
		  {
		    int c2Coarse=f2c0[pairWith];
		    if(c2Coarse==-1 && f2c0[c]==-1)
		      {
			f2c0[pairWith]=count0;
			f2c0[c]=count0;
			CoarseFineCount[count0]++;
			count0++;
		      }
		    else if(f2c0[c]==-1 && c2Coarse>=0)
		      f2c0[c]=c2Coarse;/*
		    else if(f2c0[c]!=-1 && c2Coarse==-1)
		      f2c0[pairWith]=f2c0[c];
		       */
		  }
	      }
	  }
      }
    else
      {
	for(int f=0;f<inFaceCount;f++)
	  {
	    const int c=faceCells0(f,0);
	    if(f2c0[c]<0) //dont bother if im already paired
	      {
		//loop through all neighbors to find pairing
		const int neibCount=cellCells0.getCount(c);
		pairWith=-1;
		T maxArea=0.;
		int c2;
		for(int neib=0;neib<neibCount;neib++)
		  {
		    const int fglob=globCellFaces0(c,neib);
		    
		    if(inBCfArray0[fglob]==0)  //not a boundary face
		      {
			if(c==globFaceCells0(fglob,1))
			  c2=globFaceCells0(fglob,0);
			else
			  c2=globFaceCells0(fglob,1);
			
			if(f2c0[c2]==-1)  //not already paired
			  {
			    if(areaMagArray0[fglob]>maxArea)
			      {
				pairWith=c2;
				maxArea=areaMagArray0[fglob];
			      }
			  }

		      }
		  }
		
		if(pairWith!=-1)
		  {
		    f2c0[pairWith]=count0;
		    f2c0[c]=count0;
		    count0++;
		  }
	      }
	  }
      }

    StorageSite preCoarse0(count0);
    CRPtr preCoarseToFineCells0=CRPtr(new CRConnectivity(preCoarse0,cells0));
    preCoarseToFineCells0->initCount();

    for(int c=0;c<cells0.getSelfCount();c++)
      if(f2c0[c]!=-1)
	preCoarseToFineCells0->addCount(f2c0[c],1);


    preCoarseToFineCells0->finishCount();

    for(int c=0;c<cells0.getSelfCount();c++)
      if(f2c0[c]!=-1)
	preCoarseToFineCells0->add(f2c0[c],c);

    preCoarseToFineCells0->finishAdd();

    //now use mesh0 pairing to pair mesh1
    for(int f=0;f<inFaceCount;f++)
      {
	const int f1=common01[f];
	const int c01ghost=faceCells1(f1,1);
	const int c01=faceCells1(f1,0);
	const int c00=faceCells0(f,0);
	if(f2c1[c01]<0 && f2c0[c00]>-1) //dont bother if im already paired, mesh0 must be paired
	  {
	    //const int neibCount0=cellCells0.getCount(c00);
	    const int coarseC00=f2c0[c00];
	    const int cC00fineCnt=preCoarseToFineCells0->getCount(coarseC00);
	    
	    for(int cC00fine=0;cC00fine<cC00fineCnt;cC00fine++)   //who is mesh0 paired with
	      {
		//const int c10=cellCells0(c00,neib0);
		const int fineCell=(*preCoarseToFineCells0)(coarseC00,cC00fine);
		if(fineCell!=c00)//(f2c0[c00]==f2c0[c10])
		  {
		    const int c10gNeibs=cellsIntGhost0.getCount(fineCell);  //cell is on an interface
		    for(int c10ghost=0;c10ghost<c10gNeibs;c10ghost++)
		      {
			const int c10g=cellsIntGhost0(fineCell,c10ghost);
			const int f10=cellFaces0(c10g,0);
			const int f11=common01[f10];
			const int c11ghost=faceCells1(f11,1);
			const int c11=faceCells1(f11,0);
			const int neibCount1=cellCells1.getCount(c01);
			bool isNeib(false);
			for(int c01neibs=0;c01neibs<neibCount1;c01neibs++)
			  {
			    if(cellCells1(c01,c01neibs)==c11)
			      {
				isNeib=true;
				break;
			      }
			  }

			if(isNeib)  //direct neighbor
			  {
			    if(f2c1[c11]>-1 && f2c1[c01]<0)  //other guy paired, im not
			      f2c1[c01]=f2c1[c11];
			    else if(f2c1[c01]>-1 && f2c1[c11]<0)  //im paired, other guy's not
			      f2c1[c11]=f2c1[c01];
			    else if(f2c1[c11]<0 && f2c1[c01]<0)  //both not paired
			      {
				f2c1[c11]=count1;
				f2c1[c01]=count1;
				count1++;
			      }
			  }
			else if(_level==-1)   //check for indirect neighbors
			  {
			    const int c11neibCnt=cellCells1.getCount(c11);
			    int middleMan=-1;
			    for(int c11nb=0;c11nb<c11neibCnt;c11nb++)
			      {
				const int c11neib=cellCells1(c11,c11nb);
				for(int c01neibs=0;c01neibs<neibCount1;c01neibs++)
				  {
				    if(cellCells1(c01,c01neibs)==c11neib)
				      {
					middleMan=c11neib;
					break;
				      }
				  }
			      }
			    if(middleMan!=-1)
			      {
				if(f2c1[middleMan]==-1)
				  {
				    if(f2c1[c11]<0 && f2c1[c01]<0)  //both not paired
				      {
					f2c1[c11]=count1;
					f2c1[c01]=count1;
					f2c1[middleMan]=count1;
					count1++;
				      }
				    else if(f2c1[c11]>-1 && f2c1[c01]<0)  //other guy paired, im not
				      {
					f2c1[c01]=f2c1[c11];
					f2c1[middleMan]=f2c1[c11];
				      }
				    else if(f2c1[c01]>-1 && f2c1[c11]<0)  //im paired, other guy's not
				      {
					f2c1[c11]=f2c1[c01];
					f2c1[middleMan]=f2c1[c01];
				      }
				  }
			      }
			  }

		      }
		  }
	      }

	    /*
	    if(f2c1[c01]==-1)  //still hasnt been paired -- means mesh0 cell wasnt paired.
	      {
		f2c1[c01]=count1;
		count1++;
	      }*/
	    
	  }
      }
    
  }

  void doSweeps(const int sweeps)
  {
    for(int sweepNo=0;sweepNo<sweeps;sweepNo++)
      {
	smooth(1);
	//smooth(-1);
      }
    //applyTemperatureBoundaries();
  }

  void smooth(int dir)
  {

    const int numMeshes=_meshes.size();
    int start;
    if(dir==1)
      start=0;
    else
      start=numMeshes-1;

    COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);

    for(int msh=start;((msh<numMeshes)&&(msh>-1));msh+=dir)
      {
	const Mesh& mesh=*_meshes[msh];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[msh]];
	FgTKClistMap FgToKsc;

	if(!_MeshToIC.empty())
	  {
	    IntArray& ICs=*_MeshToIC[msh];
	    const int totIC=ICs.getLength();
	    for(int i=0;i<totIC;i++)
	      {
		COMETIC<T>& ic=*_IClist[ICs[i]];
		int fgid=ic.getSelfFaceID(msh);
		FgToKsc[fgid]=ic.getKConnectivity(fgid);
	      }
	  }
	
	const BCcellArray& BCArray=*(_BCells[msh]);
	const BCfaceArray& BCfArray=*(_BFaces[msh]);

	COMETDiscretizer<T> CDisc(mesh,_geomFields,_macro,
				  kspace,_bcMap,BCArray,BCfArray,_options,
				  FgToKsc);

	CDisc.setfgFinder();

#ifdef FVM_PARALLEL
	swapGhostInfo();
#endif

	if(_options.Scattering=="Full")
	  {
	    applyTemperatureBoundaries();
	    if(_level==0)
	      {
		CDisc.COMETSolveFull(dir,_level);
		//CDisc.COMETSolveFull(-dir,_level);
	      }
	    else
	      {
		CDisc.COMETSolveCoarse(dir,_level); 
		CDisc.COMETSolveCoarse(-dir,_level); 
	      }
	  }
	else
	  {
	    if(_level>0 || _options.Convection=="FirstOrder")
	      {
		applyTemperatureBoundaries();
		CDisc.COMETSolveCoarse(dir,_level); 
		CDisc.COMETSolveCoarse(-dir,_level); 
	      }
	    else
	      {
		//applyTemperatureBoundaries();
		CDisc.COMETSolveFine(dir,_level);
		CDisc.COMETSolveFine(-dir,_level);
	      }
	  }
	
	if(!_MeshToIC.empty())
	  {
	    IntArray& ICs=*_MeshToIC[msh];
	    const int totIC=ICs.getLength();
	    for(int i=0;i<totIC;i++)
	      {
		COMETIC<T>& ic=*_IClist[ICs[i]];
		ComInt.updateOtherGhost(ic,msh,_level!=0);
	      }
	  }
	
      }

  }
  
  T updateResid(const bool addFAS)
  {
    const int numMeshes=_meshes.size();

    COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);

    const int listSize=_IClist.size();
    for(int ic=0;ic<listSize;ic++)
      {
	COMETIC<T>* icPtr=_IClist[ic];
	//if(icPtr->InterfaceModel=="DMM" && _level==0)
	//ComInt.remakeDMMcoeffs(*icPtr);
	ComInt.updateResid(*icPtr,addFAS);
      }

    T highResid=-1.;
    T currentResid;

#ifdef FVM_PARALLEL
    swapGhostInfo();
#endif

    for(int msh=0;msh<numMeshes;msh++)
      {
	const Mesh& mesh=*_meshes[msh];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[msh]];
	const BCcellArray& BCArray=*(_BCells[msh]);
	const BCfaceArray& BCfArray=*(_BFaces[msh]);
	FgTKClistMap FgToKsc;
	COMETDiscretizer<T> CDisc(mesh,_geomFields,_macro,
				  kspace,_bcMap,BCArray,BCfArray,_options,
				  FgToKsc);
	
	CDisc.setfgFinder();
	if(_options.Scattering=="Full")
	  {
	    if(_level==0)
	      CDisc.findResidFull();
	    else
	      CDisc.findResidCoarse(addFAS);
	  }
	else
	  {
	    if(_level>0 || _options.Convection=="FirstOrder")
	      CDisc.findResidCoarse(addFAS);
	    else
	      CDisc.findResidFine();
	  }

	currentResid=CDisc.getAveResid();

	if(highResid<0)
	  highResid=currentResid;
	else
	  if(currentResid>highResid)
	    highResid=currentResid;
      }

    return highResid;
  }

  void swapGhostInfo()
  {
    const int numMeshes=_meshes.size();
    for(int msh=0;msh<numMeshes;msh++)
      {
	const Mesh& mesh=*_meshes[msh];
	const StorageSite& cells=mesh.getCells();
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[msh]];
	kspace.syncLocal(cells);
      }
  }

  void cycle()
  {
    doSweeps(_options.preSweeps);
    
    if(_level+1<_options.maxLevels)
      {
	updateResid(_level!=0);
	injectResid();
	_coarserLevel->sete0();
	_coarserLevel->makeFAS();
	_coarserLevel->cycle();
	correctSolution();
      }
    
    doSweeps(_options.postSweeps);
  }

  void injectResid()
  {
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& finerMesh=*_meshes[n];
	int Knum=_MeshKspaceMap[n];
	Tkspace& finerKspace=*_kspaces[Knum];
	const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
	Tkspace& coarserKspace=_coarserLevel->getKspace(Knum);
	PhononMacro& coarserMacro=_coarserLevel->getMacro();
	const StorageSite& finerCells=finerMesh.getCells();
	const StorageSite& coarserCells=coarserMesh.getCells();
	const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);
	const TArray& coarserVol=dynamic_cast<TArray&>(_geomFields.volume[coarserCells]);
	const TArray& finerVol=dynamic_cast<TArray&>(_geomFields.volume[finerCells]);

	TArray& coarserVar=coarserKspace.geteArray();
	TArray& coarserInj=coarserKspace.getInjArray();
	TArray& coarserFAS=coarserKspace.getFASArray();
	TArray& finerVar=finerKspace.geteArray();
	TArray& finerRes=finerKspace.getResArray();

	const int cellCount=coarserCells.getSelfCount();
	const int cellTotCount=coarserCells.getCount();
	coarserVar.zero();
	coarserFAS.zero();

	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
	    
	    //sum up contributions from fine cells
	    for(int fc=0;fc<fineCount;fc++)
	      {
		const int cell=CoarserToFiner(c,fc);
		const int klen=finerKspace.getlength();
		int coarserCellIndex=coarserKspace.getGlobalIndex(c,0);
		int finerCellIndex=finerKspace.getGlobalIndex(cell,0);

		for(int k=0;k<klen;k++)
		  {
		    Tkvol& kvol=finerKspace.getkvol(k);
		    const int numModes=kvol.getmodenum();
		    for(int m=0;m<numModes;m++)
		      {
			coarserVar[coarserCellIndex]+=finerVar[finerCellIndex]*finerVol[cell];
			coarserFAS[coarserCellIndex]+=finerRes[finerCellIndex];
			finerCellIndex++;
			coarserCellIndex++;
		      }
		  }
	      }

	    //make volume average
	    int coarserCellIndex=coarserKspace.getGlobalIndex(c,0);
	    const int klen=finerKspace.getlength();
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=finerKspace.getkvol(k);
		//Tkvol& ckvol=coarserKspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    coarserVar[coarserCellIndex]/=coarserVol[c];
		    coarserInj[coarserCellIndex]=coarserVar[coarserCellIndex];
		    coarserCellIndex++;
		  }
	      }
	  }

	for(int c=cellCount;c<cellTotCount;c++)
	  {
	    const int cell=CoarserToFiner(c,0);
	    int coarserCellIndex=coarserKspace.getGlobalIndex(c,0);
	    int finerCellIndex=finerKspace.getGlobalIndex(cell,0);
	    const int klen=finerKspace.getlength();
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=finerKspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    coarserVar[coarserCellIndex]=0.;
		    coarserFAS[coarserCellIndex]=0.;
		    coarserVar[coarserCellIndex]+=finerVar[finerCellIndex];
		    coarserFAS[coarserCellIndex]+=finerRes[finerCellIndex];	
		    coarserInj[coarserCellIndex]=coarserVar[coarserCellIndex];
		    coarserCellIndex++;
		    finerCellIndex++;
		  }
	      }
	  }	
	
	TArray& coarserVarM=dynamic_cast<TArray&>(coarserMacro.temperature[coarserCells]);
	TArray& coarserInjM=dynamic_cast<TArray&>(coarserMacro.TlInjected[coarserCells]);
	TArray& coarserFASM=dynamic_cast<TArray&>(coarserMacro.TlFASCorrection[coarserCells]);
	TArray& finerVarM=dynamic_cast<TArray&>(_macro.temperature[finerCells]);
	TArray& finerResM=dynamic_cast<TArray&>(_macro.TlResidual[finerCells]);
	    
	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
	    coarserVarM[c]=0.;
	    coarserFASM[c]=0.;
		
	    for(int fc=0;fc<fineCount;fc++)
	      {
		const int cell=CoarserToFiner(c,fc);
		coarserVarM[c]+=finerVarM[cell]*finerVol[cell];
		coarserFASM[c]+=finerResM[cell];
	      }
	    coarserVarM[c]/=coarserVol[c];
	    coarserInjM[c]=coarserVarM[c];
	  }	
      }
  }
  
  void makeFAS()
  {
    updateResid(false);

    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	const StorageSite& cells=mesh.getCells();
	kspace.makeFAS();

	TArray& resArray=dynamic_cast<TArray&>(_macro.TlResidual[cells]);
	TArray& fasArray=dynamic_cast<TArray&>(_macro.TlFASCorrection[cells]);
	fasArray-=resArray;
      }
  }

  void sete0()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	const StorageSite& cells=mesh.getCells();
	const int cellCount=cells.getSelfCount();
	TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[cells]);
	TArray& e0Array=kspace.gete0Array();
	const int klen=kspace.getlength();

	for(int c=0;c<cellCount;c++)
	  {
	    const T Tlat=Tl[c];
	    int cellIndex=kspace.getGlobalIndex(c,0);
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    e0Array[cellIndex]=mode.calce0(Tlat);
		    kspace.updateTau(c,Tl[c]);
		    cellIndex++;
		  }
	      }
	  }
      }
  }

  void correctSolution()
  {
    const int numMeshes = _meshes.size();
    
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& finerMesh=*_meshes[n];
	int Knum=_MeshKspaceMap[n];
	Tkspace& finerKspace=*_kspaces[Knum];
	const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
	Tkspace& coarserKspace=_coarserLevel->getKspace(Knum);
	PhononMacro& coarserMacro=_coarserLevel->getMacro();
	const StorageSite& finerCells=finerMesh.getCells();
	const StorageSite& coarserCells=coarserMesh.getCells();
	const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);
	const CRConnectivity& coarseCellFaces=coarserMesh.getCellFaces();
	BCfaceArray& BCfArray=*(_BFaces[n]);
	
	TArray& coarserArray=coarserKspace.geteArray();
	TArray& finerArray=finerKspace.geteArray();
	TArray& injArray=coarserKspace.getInjArray();
	
	const int cellCount=coarserCells.getSelfCount();
	const int cellTotCount=coarserCells.getCount();

	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
	    for(int fc=0;fc<fineCount;fc++)
	      {
		int coarserCellIndex=coarserKspace.getGlobalIndex(c,0);
		const int finerCell=CoarserToFiner(c,fc);
		int finerCellIndex=finerKspace.getGlobalIndex(finerCell,0);
		const int klen=finerKspace.getlength();
		for(int k=0;k<klen;k++)
		  {
		    Tkvol& kvol=finerKspace.getkvol(k);
		    const int numModes=kvol.getmodenum();
		    for(int m=0;m<numModes;m++)
		      {
			const T correction=coarserArray[coarserCellIndex]-injArray[coarserCellIndex];		    
			finerArray[finerCellIndex]+=correction;
			coarserCellIndex++;
			finerCellIndex++;
		      }
		  }
	      }
	  }

	for(int c=cellCount;c<cellTotCount;c++)
	  {
	    const int f=coarseCellFaces(c,0);
	    if(BCfArray[f]==4)  //correct interface cells only
	      {
		const int fineCount=CoarserToFiner.getCount(c);
		for(int fc=0;fc<fineCount;fc++)
		  {
		    const int finerCell=CoarserToFiner(c,fc);
		    int coarserCellIndex=coarserKspace.getGlobalIndex(c,0);
		    int finerCellIndex=finerKspace.getGlobalIndex(finerCell,0);
		    const int klen=finerKspace.getlength();

		    for(int k=0;k<klen;k++)
		      {
			Tkvol& kvol=finerKspace.getkvol(k);
			const int numModes=kvol.getmodenum();
			for(int m=0;m<numModes;m++)
			  {
			    const T correction=coarserArray[coarserCellIndex]-injArray[coarserCellIndex];
			    finerArray[finerCellIndex]+=correction;
			    coarserCellIndex++;
			    finerCellIndex++;
			  }
		      }

		  }
	      }
	  }

	TArray& coarserArrayM=dynamic_cast<TArray&>(coarserMacro.temperature[coarserCells]);
	TArray& injArrayM=dynamic_cast<TArray&>(coarserMacro.TlInjected[coarserCells]);
	TArray& finerArrayM=dynamic_cast<TArray&>(_macro.temperature[finerCells]);
	
	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
	    const T correction=coarserArrayM[c]-injArrayM[c];
	    
	    for(int fc=0;fc<fineCount;fc++)
	      finerArrayM[CoarserToFiner(c,fc)]+=correction;
	  }
      }
  }
  
  void updateTL()
  {
    const int numMeshes=_meshes.size();
    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*_meshes[n];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];	
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	const T tauTot=kspace.calcTauTot();
	TArray& TL=dynamic_cast<TArray&>(_macro.temperature[cells]);
	TArray& e0Array=dynamic_cast<TArray&>(_macro.e0[cells]);
	
	for(int c=0;c<numcells;c++)
	  kspace.NewtonSolve(TL[c],e0Array[c]*tauTot);
      }
  }
  
  void advance(const int iters)
  {
    _residual=updateResid(false);
    int niters=0;
    const T absTol=_options.absTolerance;
    const int show=_options.showResidual;

    while((niters<iters) && (_residual>absTol))
      {
	cycle();
	niters++;
	_residual=updateResid(false);
	if(niters%show==0)
	  cout<<"Iteration:"<<niters<<" Residual:"<<_residual<<endl;
	
      }
    
    //applyTemperatureBoundaries();
    //calcModeTemps();
  }

  T HeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    T r(0.);
    bool found = false;
    const int n=mesh.getID();
    Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
    TArray& eArray=kspace.geteArray();
    const T DK3=kspace.getDK3();

    const T hbar=6.582119e-16;  // (eV s)

    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
	  {
	    const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
	    //const StorageSite& cells = mesh.getCells();
	    const CRConnectivity& faceCells=mesh.getFaceCells(faces);
	    const Field& areaField=_geomFields.area;
	    const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]);
	    
	    for(int f=0; f<nFaces; f++)
	      {
		const VectorT3 An=faceArea[f];
		const int c1=faceCells(f,1);
		int cellIndex=kspace.getGlobalIndex(c1,0);
		for(int k=0;k<kspace.getlength();k++)
		  {
		    Tkvol& kv=kspace.getkvol(k);
		    int modenum=kv.getmodenum();
		    for(int m=0;m<modenum;m++)
		      {
			VectorT3 vg=kv.getmode(m).getv();
			T dk3=kv.getdk3();
			T energy=hbar*kv.getmode(m).getomega();
			const T vgdotAn=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
			r += eArray[cellIndex]*vgdotAn*(dk3/DK3);//*energy;
			cellIndex++;
		      }
		  }
	      }
	    found=true;
	  }
      }

#ifdef FVM_PARALLEL
    found=true;
    int one=1;
    double tempR=r;
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempR, one, MPI::DOUBLE, MPI::SUM);
    r=tempR;
#endif

    if (!found)
      throw CException("getHeatFluxIntegral: invalid faceGroupID");
    return r*DK3;
  }

  ArrayBase* modewiseHeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    bool found = false;
    const int n=mesh.getID();
    Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
    TArray& eArray=kspace.geteArray();
    const T DK3=kspace.getDK3();
    TArray* qptr(new TArray(kspace.gettotmodes()));
    TArray& q(*qptr);
    q.zero();

    const T hbar=6.582119e-16;  // (eV s)

    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
        const FaceGroup& fg = *fgPtr;
        if (fg.id == faceGroupId)
	  {
	    const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
	    //const StorageSite& cells = mesh.getCells();
	    const CRConnectivity& faceCells=mesh.getFaceCells(faces);
	    const Field& areaField=_geomFields.area;
	    const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]);
	    
	    for(int f=0; f<nFaces; f++)
	      {
		const VectorT3 An=faceArea[f];
		const int c1=faceCells(f,1);
		int cellIndex=kspace.getGlobalIndex(c1,0);
		for(int k=0;k<kspace.getlength();k++)
		  {
		    Tkvol& kv=kspace.getkvol(k);
		    int modenum=kv.getmodenum();
		    for(int m=0;m<modenum;m++)
		      {
			VectorT3 vg=kv.getmode(m).getv();
			const int index=kv.getmode(m).getIndex()-1;
			T dk3=kv.getdk3();
			T energy=hbar*kv.getmode(m).getomega();
			const T vgdotAn=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
			q[index]+= eArray[cellIndex]*vgdotAn*dk3;//*energy;
			cellIndex++;
		      }
		  }
	      }
	    found=true;
	  }
      }
    
    return qptr;
  }

  T getWallArea(const Mesh& mesh, const int faceGroupId)
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
	    const Field& areaMagField=_geomFields.areaMag;
	    const TArray& faceArea=dynamic_cast<const TArray&>(areaMagField[faces]);
	    	    
       	    for(int f=0; f<nFaces; f++)
	      r += faceArea[f];
	    
	    found=true;
	    break;
	  }
      }

    if (!found)
      throw CException("getwallArea: invalid faceGroupID");
    return r;
  }

  VectorT3 getWallAreaVector(const Mesh& mesh, const int faceGroupId)
  {
    VectorT3 An;
    An=0.;
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;
	if (fg.id == faceGroupId)
	  {
	    const StorageSite& faces = fg.site;
	    const int nFaces = faces.getCount();
	    const Field& areaField=_geomFields.area;
	    const VectorT3Array& faceArea=dynamic_cast<const VectorT3Array&>(areaField[faces]);
	    for(int f=0; f<nFaces; f++)
	      An+=faceArea[f];
	    found=true;
	    break;
	  }
      }

#ifdef FVM_PARALLEL
    found=true;
    int one=1;
    double tempR=An[0];
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempR, one, MPI::DOUBLE, MPI::SUM);
    An[0]=tempR;
    tempR=An[1];
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempR, one, MPI::DOUBLE, MPI::SUM);
    An[1]=tempR;
    tempR=An[2];
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempR, one, MPI::DOUBLE, MPI::SUM);
    An[2]=tempR;
#endif

    if (!found)
      throw CException("getwallArea: invalid faceGroupID");
    return An;
  }
  
  ArrayBase* getValueArray(const Mesh& mesh, const int cell)
    {
      //only returns the e" values, not the lattice temperature
      const int n=mesh.getID();
      Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
      const int allModes=kspace.gettotmodes();
      TArray* vals=new TArray(allModes);
      const StorageSite& cells=mesh.getCells();
      const int len=kspace.getlength();
      int count=0;
      for(int k=0;k<len;k++)
	{
	  Tkvol& kvol=kspace.getkvol(k);
	  const int modes=kvol.getmodenum();
	  for(int m=0;m<modes;m++)
	    {
	      Tmode& mode=kvol.getmode(m);
	      Field& eField=mode.getfield();
	      const TArray& eArray=dynamic_cast<const TArray&>(eField[cells]);
	      (*vals)[count]=eArray[cell];
	      count++;
	    }
	}
      return vals;
    }

  void equilibrate()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	const StorageSite& cells=mesh.getCells();
	const int numcells=cells.getCount();
	const TArray& Tl=dynamic_cast<const TArray&>(_macro.temperature[cells]);
	const int len=kspace.getlength();
	
	for(int k=0;k<len;k++)
	  {
	    Tkvol& kvol=kspace.getkvol(k);
	    const int modes=kvol.getmodenum();
	    for(int m=0;m<modes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Field& eField=mode.getfield();
		TArray& eArray=dynamic_cast<TArray&>(eField[cells]);
		for(int c=0;c<numcells;c++)
		  eArray[c]=mode.calce0(Tl[c]);
	      }
	  }
      }
  }

  ArrayBase& getLatticeTemp(const Mesh& mesh)
  {
    const StorageSite& cells=mesh.getCells();
    TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[cells]);
    return Tl;
  }

  void setStraightLine(const T T1, const T T2)
  {
    const int numMeshes = _meshes.size();

    T xmax(0.);
    T xmin(0.);

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	const int numcells=cells.getCount();
	const VectorT3Array& coords=dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);

	for(int c=0;c<numcells;c++)
	  {
	    T x=coords[c][0];
	    if(x>xmax)
	      xmax=x;
	    if(x<xmin)
	      xmin=x;
	  }

      }
    
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	//const int numcells=cells.getCount();
	const int selfcells=cells.getSelfCount();
	TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[cells]);
	const VectorT3Array& coords=dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
	
	for(int c=0;c<selfcells;c++)
	  {
	    T factor=(coords[c][0]-xmin)/(xmax-xmin);
	    Tl[c]=T1+factor*(T2-T1);
	  }

      }

  }

  void calcModeTemps()
  {
    const int numMeshes=_meshes.size();
    for (int n=0;n<numMeshes;n++)
       {
	 Mesh& mesh=*_meshes[n];
	 if(_MeshKspaceMap[n]==-1)
	   throw CException("Have not set the Kspace for this Mesh!!");
	 Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	 const int numK=kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 FieldVector& FieldVec=*_macro.BranchTemperatures[n];
	 const int modeCount=FieldVec.size();
	 TArray eSum(numcells);
	 const TArray& TL=dynamic_cast<const TArray&>(_macro.temperature[cells]);
	 
	 for(int m=0;m<modeCount;m++)
	   {
	     eSum.zero();
	     Field& modeField=*FieldVec[m];
	     TArray& modeTemp=dynamic_cast<TArray&>(modeField[cells]);

	     for(int k=0;k<numK;k++)
	       {
		 T dk3=kspace.getkvol(k).getdk3();
		 Tmode& mode=kspace.getkvol(k).getmode(m);
		 Field& eField=mode.getfield();
		 TArray& eArray=dynamic_cast<TArray&>(eField[cells]);
		 
		 for(int c=0;c<numcells;c++)
		   eSum[c]+=eArray[c]*dk3;
	       }

	     for(int c=0;c<numcells;c++)
	       modeTemp[c]=kspace.calcModeTemp(TL[c],eSum[c],m);
	   }

       }
  }

  void calcModeFlux()
  {
    const int numMeshes=_meshes.size();
    const T eVtoJoule=1.60217646e-19;
    for (int n=0;n<numMeshes;n++)
       {
	 Mesh& mesh=*_meshes[n];
	 if(_MeshKspaceMap[n]==-1)
	   throw CException("Have not set the Kspace for this Mesh!!");
	 Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	 const int numK=kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 const int modeCount=kspace.getkvol(0).getmodenum();
	 VectorT3Array& Q=dynamic_cast<VectorT3Array&>(_macro.heatFlux[cells]);

	 FieldVector* FieldVecPtr=new FieldVector();
	 for(int m=0;m<modeCount;m++)
	   {
	     shared_ptr<Field> modeField(new Field("mode"));
	     VT3Ptr qptr(new VectorT3Array(numcells));
	     VectorT3Array& q(*qptr);
	     q.zero();
	     modeField->addArray(cells,qptr);
	     FieldVecPtr->push_back(modeField);
	     _macro.BranchFlux[n]=FieldVecPtr;

	     for(int c=0;c<numcells;c++)
	       {
		 for(int k=0;k<numK;k++)
		   {
		     Tmode& mode=kspace.getkvol(k).getmode(m);
		     T dk3=kspace.getkvol(k).getdk3();
		     const int index=mode.getIndex()-1;
		     q[c]+=mode.getv()*dk3*kspace.gete(c,index)*eVtoJoule;
		     Q[c]+=mode.getv()*dk3*kspace.gete(c,index)*eVtoJoule;
		   }

	       }
	   }

       }
  }

  T calcDomainStats()
  {
    const int numMeshes=_meshes.size();
    TArray meshVols(numMeshes);
    T totalVol(0.);
    meshVols.zero();
    T interfaceArea(0.);
    T interfaceAreaX(0.);
    T interfaceAreaY(0.);
    T interfaceAreaZ(0.);
    int interfaceFaces(0);

    for (int n=0;n<numMeshes;n++)
      {
	 Mesh& mesh=*_meshes[n];
	 const StorageSite& cells=mesh.getCells();
	 const  TArray& vol=dynamic_cast<const TArray&>(_geomFields.volume[cells]);
	 
	 for(int c=0;c<cells.getSelfCount();c++)
	   meshVols[n]+=vol[c];

	 totalVol+=meshVols[n];

	 foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	   {
	     const FaceGroup& fg = *fgPtr;
	     if(fg.id>0)
	       {
		 if(_bcMap[fg.id]->bcType == "Interface")
		   {
		     const StorageSite& faces=fg.site;
		     const int numFaces=faces.getCount();
		     const Field& AreaMagField=_geomFields.areaMag;
		     const TArray& AreaMag=
		       dynamic_cast<const TArray&>(AreaMagField[faces]);
		     const VectorT3Array& AreaDir=
		       dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
		     
		     interfaceFaces+=numFaces;

		     for(int f=0;f<numFaces;f++)
		       {
			 interfaceArea+=AreaMag[f];
			 interfaceAreaX+=fabs(AreaDir[f][0]);
			 interfaceAreaY+=fabs(AreaDir[f][1]);
			 interfaceAreaZ+=fabs(AreaDir[f][2]);
		       }

		   }
	       }
	   }
      }

    interfaceArea/=2.;
    interfaceAreaX/=2.;
    interfaceAreaY/=2.;
    interfaceAreaZ/=2.;
    interfaceFaces/=2;

    for (int n=0;n<numMeshes;n++)
      {
	Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	cout<<endl;
	cout<<"Mesh: "<<n<<endl;
	cout<<"Self Cell Count: "<<cells.getSelfCount()<<endl;
	cout<<"Ghost Cell Count: "<<cells.getCount()-cells.getSelfCount()<<endl;
	cout<<"Mesh Volume: "<<meshVols[n]<<endl;
      }

    cout<<endl;
    cout<<"Total Volume: "<<totalVol<<endl;
    cout<<"Interface Area: "<<interfaceArea<<endl;
    cout<<"Interface Area (X): "<<interfaceAreaX<<endl;
    cout<<"Interface Area (Y): "<<interfaceAreaY<<endl;
    cout<<"Interface Area (Z): "<<interfaceAreaZ<<endl;
    cout<<"Area/Volume: "<<interfaceArea/totalVol<<endl;
    cout<<"Volume/Area: "<<totalVol/interfaceArea<<endl;
    cout<<"Interface Face Count: "<<interfaceFaces<<endl;
    return interfaceArea/totalVol;
    
  }

  void makeCellColors(const int level)
  {
    
    TCOMET* thismodel=getModelPointer(level);
    MeshList& meshes=const_cast<MeshList&>(thismodel->getMeshList());
    const int numMeshes = meshes.size();
    shared_ptr<Field> colorFieldPtr(new Field("Colors"));
    _macro.CellColors.push_back(colorFieldPtr);
    Field& colorField=*colorFieldPtr;
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*meshes[n];
	const StorageSite& cells=mesh.getCells();
	const int cellCount=cells.getSelfCount();
	const int cellTotCount=cells.getCount();
	const CRConnectivity& cellCells=mesh.getCellCells();
	
	TArrptr colorPtr(new TArray(cellTotCount));
	colorField.addArray(cells, colorPtr);
	TArray& colorArray=(*colorPtr);
	colorArray=-1;

	for(int c=0;c<cellCount;c++)
	  {
	    T color=0.;
	    const int neibs=cellCells.getCount(c);

	    while(colorArray[c]==-1)
	      {

		bool same(false);

		for(int j=0;j<neibs;j++)
		  {
		    if(color==colorArray[cellCells(c,j)])
		      {
			same=true;
			break;
		      }
		  }

		if(same)
		  color++;
		else
		  colorArray[c]=color;

	      }

	  }
      }

    //if(_coarserLevel!=NULL)
    //_coarserLevel->makeCellColors();

  }

  void makePlotColors(const int level)
  {

    if(level==0)
      {
	if(_options.maxLevels>1)
	  makeFinestToCoarseConn();
	shared_ptr<Field> field0ptr(new Field("plotColor"));
	_macro.plottingCellColors.push_back(field0ptr);
	Field& field0=*field0ptr;
	Field& colorField=_macro.getColorField(0);

	const int numMeshes = _meshes.size();
	for (int n=0; n<numMeshes; n++)
	  {
	    const Mesh& mesh=*_meshes[n];
	    const StorageSite& cells=mesh.getCells();
	    //const int cellCount=cells.getSelfCount();
	    const int cellTotCount=cells.getCount();
	    //const CRConnectivity& cellCells=mesh.getCellCells();
	    TArray& colorArray=dynamic_cast<TArray&>(colorField[cells]);

	    TArrptr colorPtr(new TArray(cellTotCount));
	    (*colorPtr)=colorArray;
	    field0.addArray(cells, colorPtr);
	  }
	//if(_coarserLevel!=NULL)
	//  _coarserLevel->makePlotColors();
      }
    else
      {
	shared_ptr<Field> fieldptr(new Field("plotColor"));
	_macro.plottingCellColors.push_back(fieldptr);
	Field& field=*fieldptr;

	Field& colorField=_macro.getColorField(level);

	TCOMET* coarseModel=getModelPointer(level);
	const MeshList& coarseMeshes=coarseModel->getMeshList();
	const int numMeshes = _meshes.size();
	for(int n=0;n<numMeshes;n++)
	  {
	    const Mesh& finestMesh=*_meshes[n];
	    const Mesh& coarseMesh=*coarseMeshes[n];
	    const StorageSite& coarseCells=coarseMesh.getCells();
	    const StorageSite& finestCells=finestMesh.getCells();
	    const int finestCellTotCount=finestCells.getCount();
	    const int finestCellCount=finestCells.getSelfCount();
	    const CRConnectivity& finestToCoarse=finestMesh.getConnectivity(finestCells,coarseCells);
	    	    
	    TArrptr plotColorPtr(new TArray(finestCellTotCount));
	    TArray& plotColorArray=*plotColorPtr;
	    field.addArray(finestCells, plotColorPtr);

	    TArray& colorArray=dynamic_cast<TArray&>(colorField[coarseCells]);

	    for(int c=0;c<finestCellCount;c++)
	      plotColorArray[c]=colorArray[finestToCoarse(c,0)];
	  }

	//if(_coarserLevel!=NULL)
	//  _coarserLevel->makePlotColors();
      }
  }

  void makeFinestToCoarseConn()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	
	Mesh& mesh1=*(_coarserLevel->getMeshList()[n]);
	const StorageSite& cells1=mesh1.getCells();

	const CRConnectivity& firstCRC=mesh1.getConnectivity(cells1,cells);
	CRPtr FineToCoarseOld=firstCRC.getTranspose();

	mesh.setConnectivity(cells,cells1,FineToCoarseOld);

	for(int lvl=2;lvl<_options.maxLevels;lvl++)
	  {
	    TCOMET* coarseModel=getModelPointer(lvl);
	    Mesh& coarseMesh=*(coarseModel->getMeshList()[n]);
	    const StorageSite& coarseCells=coarseMesh.getCells();
	    TCOMET* finerModel=coarseModel->getFinerModel();
	    Mesh& finerMesh=*(finerModel->getMeshList()[n]);
	    const StorageSite& finerCells=finerMesh.getCells();
	    const CRConnectivity& coarseToFine=coarseMesh.getConnectivity(coarseCells, finerCells);
	    CRPtr FineToCoarseNew=coarseToFine.getTranspose();
	    FineToCoarseNew=FineToCoarseOld->multiply(*FineToCoarseNew,true);
	    mesh.setConnectivity(cells, coarseCells,FineToCoarseNew);
	    FineToCoarseOld=FineToCoarseNew;
	  }
	
      }
  }

  TCOMET* getModelPointer(const int level)
  {
    if(_level==level)
      return this;
    else
      return _coarserLevel->getModelPointer(level);
  }

  bool sameFaceGroup(const Mesh& mesh, const int f0, const int f1)
  {
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg=*fgPtr;
	const int off=fg.site.getOffset();
	const int cnt=fg.site.getCount();
	const int fin=off+cnt;

	if(f0<fin && f0>=off)
	  if(f1<fin && f1>=off)
	    return true;
      }
    return false;
  }

  T getAverageTemperature()
  {
    const int numMeshes = _meshes.size();
    T r(0);
    T volTot(0);
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	const int numcells=cells.getSelfCount();
	const TArray& Tl=dynamic_cast<const TArray&>(_macro.temperature[cells]);
	const TArray& cv=dynamic_cast<const TArray&>(_geomFields.volume[cells]);
	for(int c=0;c<numcells;c++)
	  {
	    r+=Tl[c]*cv[c];
	    volTot+=cv[c];
	  }
      }
    return r/volTot;
  }

  void makeNonEqTemp()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	const int numcells=cells.getSelfCount();
	TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[cells]);
	Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
	const TArray& e=kspace.geteArray();

	const int len=kspace.getlength();

	for(int c=0;c<numcells;c++)
	  Tl[c]=kspace.calcLatTemp(c);	    

      }
  }
  
  void setBCMap(COMETBCMap* bcMap) {_bcMap=*bcMap;}
  void setCoarserLevel(TCOMET* cl) {_coarserLevel=cl;}
  void setFinerLevel(TCOMET* fl) {_finerLevel=fl;}
  TCOMET* getFinerModel() {return _finerLevel;}
  int getLevel() {return _level;}
  const MeshList& getMeshList() {return _meshes;}
  GeomFields& getGeomFields() {return _geomFields;}
  TkspList& getKspaces() {return _kspaces;}
  Tkspace& getKspace(const int i) {return *_kspaces[i];}
  PhononMacro& getMacro() {return _macro;}
  T getResidual() 
  {
#ifdef FVM_PARALLEL
    int one=1;
    double tempResid=_residual;
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &tempResid, one, MPI::DOUBLE, MPI::SUM);
    _residual=tempResid;
#endif
    return _residual;
  }
  IClist& getIClist() {return _IClist;}
  MeshICmap& getMeshICmap() {return _MeshToIC;}

 private:

  const int _level;    //0 being the finest level
  GeomFields& _geomFields;
  TkspList& _kspaces;
  PhononMacro& _macro;
  TCOMET* _finestLevel;
  TCOMET* _coarserLevel;
  TCOMET* _finerLevel;
  COMETBCMap _bcMap;
  TCModOpts _options;
  BCfaceList _BFaces;
  BCcellList _BCells;
  T _residual;
  MeshKspaceMap _MeshKspaceMap;
  IClist _IClist;
  MeshICmap _MeshToIC;
  int _rank;

};

#endif
