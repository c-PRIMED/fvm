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
  
 COMETModel(const MeshList& meshes, const int level, GeomFields& geomFields,
	    TkspList& kspaces, PhononMacro& macro):
  Model(meshes),
    _level(level),
    _geomFields(geomFields),
    _kspaces(kspaces),
    _macro(macro),
    _residual(1.0),
    _MeshToIC()
      {
	const int numMeshes = _meshes.size();
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
		foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
		  {
		    const FaceGroup& fg = *fgPtr;
		    if (_bcMap.find(fg.id) == _bcMap.end() && fg.id!=0)
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
	 BCfaceArray& BCfArray=*(_BFaces[n]);
	 BCcellArray& BCArray=*(_BCells[n]);
	 const int numK=kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 VectorT3 lamTemp;
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 shared_ptr<TArray> deltaTcell(new TArray(numcells));
	 shared_ptr<VectorT3Array> lamArray(new VectorT3Array(numcells));
	 lamTemp.zero();
	 *lamArray=lamTemp;
	 *deltaTcell=0.;
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 _macro.deltaT.addArray(cells,deltaTcell);
	 _macro.lam.addArray(cells,lamArray);

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
	 
	 for (int k=0;k<numK;k++)
	   {
	     Tkvol& kv=kspace.getkvol(k);
	     const int numM=kv.getmodenum();
	     const T dk3=kv.getdk3();

	     for (int m=0;m<numM;m++)
	       {
		 Tmode& mode=kv.getmode(m);
		 Field& efield=mode.getfield();
		 Field& e0field=mode.gete0field();
		 Field& resfield=mode.getresid();
		 TArrptr e0var=shared_ptr<TArray>(new TArray(numcells));
		 TArrptr evar=shared_ptr<TArray>(new TArray(numcells));
		 TArrptr resid=shared_ptr<TArray>(new TArray(numcells));
		 const T einit=mode.calce0(Tinit);
		 *evar=einit;
		 *e0var=einit;
		 *resid=0.;
		 efield.addArray(cells,evar);
		 e0field.addArray(cells,e0var);
		 resfield.addArray(cells,resid);

		 if(_options.withNormal)
		   {
		     TArrptr eShifted=shared_ptr<TArray>(new TArray(numcells));
		     *eShifted=einit;
		     Field& eShiftedField=mode.geteShifted();
		     eShiftedField.addArray(cells,eShifted);
		   }
		 
		 Refl_Map& rmap=mode.getreflmap();
		 VectorT3 vg=mode.getv();
		 T vmag=sqrt(pow(vg[0],2)+pow(vg[1],2)+pow(vg[2],2));
		 VectorT3 si=vg/vmag;
		 VectorT3 so;
		 
		 //reflections
		 foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
		   {
		     const FaceGroup& fg = *fgPtr;
		     if(fg.id!=0)
		       {
			 if(_bcMap[fg.id]->bcType == "reflecting")
			   {
			     const StorageSite& faces = fg.site;
			     const Field& AreaMagField=_geomFields.areaMag;
			     const TArray& AreaMag=
			       dynamic_cast<const TArray&>(AreaMagField[faces]);
			     const Field& AreaDirField=_geomFields.area;
			     const VectorT3Array& AreaDir=
			       dynamic_cast<const VectorT3Array&>(AreaDirField[faces]);
			     
			     const VectorT3 n=AreaDir[0]/AreaMag[0];
			     const T sidotn=si[0]*n[0]+si[1]*n[1]+si[2]*n[2];
			     
			     if (sidotn > T_Scalar(0.0))
			       {
				 so=si-2.*(si[0]*n[0]+si[1]*n[1]+si[2]*n[2])*n;
				 T soMag=sqrt(pow(so[0],2)+pow(so[1],2)+pow(so[2],2));
				 so/=soMag;
				 so*=vmag;
				 Refl_pair refls;
				 Refl_pair reflsFrom;
				 kspace.findSpecs(dk3,vmag,m,so,refls);
				 rmap[fg.id]=refls;
				 const int k1=refls.first.second;
				 Tmode& mode2=kspace.getkvol(k1).getmode(m);
				 Refl_Map& rmap2=mode2.getreflmap();
				 reflsFrom.first.second=-1;
				 reflsFrom.second.second=k;
				 rmap2[fg.id]=reflsFrom;
			       }
			     else
			       {
				 Refl_pair refls;
				 refls.first.second=-1;
				 refls.second.second=-1;
				 rmap[fg.id]=refls;
			       }
			   }
		       }
		   }
	       }
	   }
	 
	 
	 shared_ptr<TArray> TlResidCell(new TArray(numcells));
	 *TlResidCell=0.;
	 _macro.TlResidual.addArray(cells,TlResidCell);

	 //setting facegroups
	 foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     if(fg.id!=0)
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
	 /*
	 const int mid0=icPtr->MeshID0;
	 const int mid1=icPtr->MeshID1;
	 ComInt.updateOtherGhost(*icPtr,mid0,false);
	 ComInt.updateOtherGhost(*icPtr,mid1,false);
	 */
       }
     
     applyTemperatureBoundaries();
     _residual=updateResid(false);
     cout<<"Creating Coarse Levels..."<<endl;
     MakeCoarseModel(this);
     cout<<"Coarse Levels Completed."<<endl;
  }

  void initFromOld()
  {
    const int numMeshes=_meshes.size();

     for (int n=0;n<numMeshes;n++)
       {
	 Mesh& mesh=*_meshes[n];

	 BCfaceArray& BCfArray=*(_BFaces[n]);
	 BCcellArray& BCArray=*(_BCells[n]);
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();

	 //setting facegroups
	 foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     if(fg.id!=0)
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
	 /*
	 const int mid0=icPtr->MeshID0;
	 const int mid1=icPtr->MeshID1;
	 ComInt.updateOtherGhost(*icPtr,mid0,false);
	 ComInt.updateOtherGhost(*icPtr,mid1,false);
	 */
       }
     
     applyTemperatureBoundaries();
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
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 shared_ptr<TArray> deltaTcell(new TArray(numcells));
	 *deltaTcell=0.; 
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 _macro.deltaT.addArray(cells,deltaTcell);
	 
	 for (int k=0;k<numK;k++)
	   {
	     Tkvol& kv=kspace.getkvol(k);
	     const int numM=kv.getmodenum();

	     for (int m=0;m<numM;m++)
	       {
		 Tmode& mode=kv.getmode(m);
		 Field& efield=mode.getfield();
		 Field& eInjField=mode.getInjected();
		 Field& e0field=mode.gete0field();
		 Field& resfield=mode.getresid();
		 Field& FASfield=mode.getFASfield();
		 TArrptr e0var=TArrptr(new TArray(numcells));
		 TArrptr evar=TArrptr(new TArray(numcells));
		 TArrptr resid=TArrptr(new TArray(numcells));
		 TArrptr FAS=TArrptr(new TArray(numcells));
		 TArrptr inj=TArrptr(new TArray(numcells));
		 const T einit=mode.calce0(Tinit);
		 *evar=einit;
		 *inj=einit;
		 *e0var=einit;
		 *resid=0.;
		 *FAS=0;
		 efield.addArray(cells,evar);
		 eInjField.addArray(cells,inj);
		 e0field.addArray(cells,e0var);
		 resfield.addArray(cells,resid);
		 FASfield.addArray(cells,FAS);
		 
	       }
	   }

	 BCfaceArray& BCfArray=*(_BFaces[n]);
	 BCcellArray& BCArray=*(_BCells[n]);
	 foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	   {
	     FaceGroup& fg = *fgPtr;
	     if(fg.id!=0)
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
	 //const int mid0=icPtr->MeshID0;
	 //const int mid1=icPtr->MeshID1;
	 //ComInt.updateOtherGhost(*icPtr,mid0);
	 //ComInt.updateOtherGhost(*icPtr,mid1);
       }
     
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
	(*scatter10)[i]=c10;
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
	foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    if(fg.id!=0)
	      {   
		const StorageSite& faces = fg.site;
		const COMETBC<T>& bc = *_bcMap[fg.id];
		
		COMETBoundary<T> cbc(faces, mesh,_geomFields,kspace,_options,fg.id);
		
		if(bc.bcType=="temperature")
		  {	      
		    FloatValEvaluator<T>
		      bTemperature(bc.getVal("specifiedTemperature"),faces);
		    
		    cbc.applyTemperatureWall(bTemperature);
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

	    int newCount= MakeCoarseMesh(finerModel->getMeshList(),
					 finerModel->getGeomFields(),
					 *newMeshesPtr);
	    	    
	    TCOMET* newModelPtr=new COMETModel(*newMeshesPtr,thisLevel,
					       finerModel->getGeomFields(),
					       *newKspacesPtr,*newMacroPtr);
	    
	    newModelPtr->setFinerLevel(finerModel);
	    finerModel->setCoarserLevel(newModelPtr);
	    newModelPtr->getOptions()=finerModel->getOptions();
	    newModelPtr->getBCs()=finerModel->getBCs();
	    newModelPtr->getMKMap()=finerModel->getMKMap();

	    newModelPtr->initCoarse();
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

  int MakeCoarseMesh(const MeshList& inMeshes, GeomFields& inGeomFields,
		      MeshList& outMeshes)
  {

    int smallestMesh=-1;
    const int numMeshes=inMeshes.size();
    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*inMeshes[n];
	const int dim=mesh.getDimension();
	Mesh* newMeshPtr=new Mesh(dim);
	outMeshes.push_back(newMeshPtr);
	newMeshPtr->setID(n);

	const StorageSite& inCells=mesh.getCells();
	StorageSite& outCells=newMeshPtr->getCells();
	StorageSite& outFaces=newMeshPtr->getFaces();
	const StorageSite& inFaces=mesh.getFaces();
	const int inCellCount=inCells.getSelfCount();
	const int inCellTotal=inCells.getCount();
	const int inFaceCount=inFaces.getCount();
	const int inGhost=inCellTotal-inCellCount;
	int coarseCount=0;
	IntArray FineToCoarse(inCellTotal);
	FineToCoarse=-1;

	const CRConnectivity& inCellinFaces=mesh.getCellFaces();
	const CRConnectivity& inFaceinCells=mesh.getFaceCells(inFaces);
	Field& areaMagField=inGeomFields.areaMag;
	const TArray& areaMagArray=dynamic_cast<const TArray&>(areaMagField[inFaces]);
	const BCfaceArray& inBCfArray=*(_BFaces[n]);

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
			
			if(FineToCoarse[c2]==-1)
			  if(areaMagArray[f]>maxArea)
			    pairWith=c2;
		      }
		  }
		
		if(pairWith!=-1)
		  {
		    FineToCoarse[c]=coarseCount;
		    FineToCoarse[c2]=coarseCount;
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
		pairWith=-1;
		
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
			  }
		      }
		  }
		
		if(pairWith==-1)
		  {
		    FineToCoarse[c]=coarseCount;
		    coarseCount++;
		  }
		else
		  {
		    if(FineToCoarse[c2perm]==-1)
		      {
			FineToCoarse[c]=coarseCount;
			FineToCoarse[c2perm]=coarseCount;
			coarseCount++;
		      }
		    else
		      FineToCoarse[c]=pairWith;
		  }
	      }	    
	  }
	
	int coarseGhost=coarseCount;
	for(int c=inCellCount;c<inCellTotal;c++)
	  {
	    FineToCoarse[c]=coarseGhost;
	    coarseGhost++;
	  }
	
	//make the coarse cell to fine cell connectivity.
	outCells.setCount(coarseCount,inGhost);
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
	int fc0,fc1,cc0,cc1;
	for(int f=0;f<inFaceCount;f++)
	  {
	    fc0=inFaceinCells(f,0);
	    fc1=inFaceinCells(f,1);
	    cc0=FineToCoarse[fc0];
	    cc1=FineToCoarse[fc1];
	    if(cc0!=cc1)
	      {
		FineFacesCoarseCells->add(f,cc0);
		FineFacesCoarseCells->add(f,cc1);
	      }
	  }

	FineFacesCoarseCells->finishAdd();

	CRPtr CoarseCellsFineFaces=FineFacesCoarseCells->getTranspose();
	CRPtr CellCellCoarse=CoarseCellsFineFaces->multiply(*FineFacesCoarseCells,true);

	int counter=0;
	BArray counted(outCells.getCount());
	counted=false;
	for(int c=0;c<outCells.getCount();c++)
	  {
	    counted[c]=true;
	    const int neibs=CellCellCoarse->getCount(c);
	    for(int n=0;n<neibs;n++)
	      {
		const int c1=(*CellCellCoarse)(c,n);
		if(!counted[c1])
		  counter++;
	      }
	  }

	outFaces.setCount(counter);

	CRPtr CoarseCellCoarseFace=CRPtr(new CRConnectivity(outCells,outFaces));
	CoarseCellCoarseFace->initCount();

	for(int c=0;c<outCells.getCount();c++)
	  {
	    const int neibs=CellCellCoarse->getCount(c);
	    CoarseCellCoarseFace->addCount(c,neibs);
	  }

	CoarseCellCoarseFace->finishCount();

	//make cell connectivity to interior faces.
	IntArray neibCounter(outCells.getCount());
	neibCounter=0;
	counter=0;
	counted=false;
	for(int c=0;c<outCells.getSelfCount();c++)
	  {
	    counted[c]=true;
	    const int neibs=CellCellCoarse->getCount(c);
	    for(int n=0;n<neibs;n++)
	      {
		const int c1=(*CellCellCoarse)(c,n);
		if(!counted[c1] && c1<outCells.getSelfCount())
		  {
		    CoarseCellCoarseFace->add(c,counter);
		    CoarseCellCoarseFace->add(c1,counter);
		    counter++;
		    neibCounter[c]++;
		    neibCounter[c1]++;
		  }
	      }
	  }

	//make cell connectivity to boundary faces.
	for(int c=outCells.getSelfCount();c<outCells.getCount();c++)
	  {
	    const int c1=(*CellCellCoarse)(c,0);
	    CoarseCellCoarseFace->add(c1,counter);
	    CoarseCellCoarseFace->add(c,counter);
	    counter++;
	    neibCounter[c]++;
	    neibCounter[c1]++;
	  }

	CoarseCellCoarseFace->finishAdd();
	
	CRPtr CoarseFaceCoarseCell=CoarseCellCoarseFace->getTranspose();

	newMeshPtr->setConnectivity(outCells,outFaces,CoarseCellCoarseFace);
	newMeshPtr->setConnectivity(outFaces,outCells,CoarseFaceCoarseCell);

	CRPtr CoarseFacesFineFaces=CRPtr(new CRConnectivity(outFaces,inFaces));
	CoarseFacesFineFaces->initCount();

	for(int f=0;f<inFaceCount;f++)
	  {
	    int fc0=inFaceinCells(f,0);
	    int fc1=inFaceinCells(f,1);
	    const int cc0=FineToCoarse[fc0];
	    const int cc1=FineToCoarse[fc1];

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
			CoarseFacesFineFaces->addCount(face,1);
			break;
		      }
		  }
	      }
	  }

	CoarseFacesFineFaces->finishCount();

	for(int f=0;f<inFaceCount;f++)
	  {
	    int fc0=inFaceinCells(f,0);
	    int fc1=inFaceinCells(f,1);
	    const int cc0=FineToCoarse[fc0];
	    const int cc1=FineToCoarse[fc1];
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
			CoarseFacesFineFaces->add(face,f);
			break;
		      }
		  }
	      }
	  }

	CoarseFacesFineFaces->finishAdd();

	const int interiorCount=outFaces.getCount()-inGhost;
	
	const StorageSite& interiorFaces=newMeshPtr->createInteriorFaceGroup(interiorCount);

	int inOffset=interiorCount;
	foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
	  {
	    const FaceGroup& fg=*fgPtr;
	    if(fg.id!=0)
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

	//now make the geom fields
	const int outCellsCount=outCells.getSelfCount();
	VT3Ptr outCellCoordsPtr=VT3Ptr(new VectorT3Array(outCellsCount));
	TArrptr outCellVolumePtr=TArrptr(new TArray(outCellsCount));
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
	
	VolumeField.addArray(outCells,outCellVolumePtr);
	inGeomFields.coordinate.addArray(outCells,outCellCoordsPtr);

	const int outFacesCount=outFaces.getCount();
	VT3Ptr outFaceAreaPtr=VT3Ptr(new VectorT3Array(outFacesCount));
	VectorT3Array& outFA=*outFaceAreaPtr;
	TArrptr outFaceAreaMagPtr=TArrptr(new TArray(outFacesCount));
	TArray& outFAMag=*outFaceAreaMagPtr;

	Field& FaceAreaField=inGeomFields.area;
	const VectorT3Array& inFA=
	  dynamic_cast<const VectorT3Array&>(FaceAreaField[inFaces]);

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
	cout<<"Level: "<<_level+1<<" Mesh: "<<n<<" Cells: "<<outCells.getSelfCount()<<endl;

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

  void doSweeps(const int sweeps)
  {
    for(int sweepNo=0;sweepNo<sweeps;sweepNo++)
      {
	smooth(1);
	smooth(-1);
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

	CDisc.COMETSolve(dir,_level); //forward
      }

    COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);
    
    const int listSize=_IClist.size();
    for(int ic=0;ic<listSize;ic++)
      {
	COMETIC<T>* icPtr=_IClist[ic];
	const int mid0=icPtr->MeshID0;
	const int mid1=icPtr->MeshID1;
	if(icPtr->InterfaceModel=="DMM")
	  ComInt.makeDMMcoeffs(*icPtr);
	ComInt.updateOtherGhost(*icPtr,mid0,_level!=0);
	ComInt.updateOtherGhost(*icPtr,mid1,_level!=0);
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
	//const int mid0=icPtr->MeshID0;
	//const int mid1=icPtr->MeshID1;
	if(icPtr->InterfaceModel=="DMM")
	  ComInt.makeDMMcoeffs(*icPtr);
	ComInt.updateResid(*icPtr,addFAS);
	//ComInt.updateOtherGhost(*icPtr,mid0);
	//ComInt.updateOtherGhost(*icPtr,mid1);
      }
    
    T highResid=-1.;
    T currentResid;
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
	CDisc.findResid(addFAS);
	currentResid=CDisc.getAveResid();

	if(highResid<0)
	  highResid=currentResid;
	else
	  if(currentResid>highResid)
	    highResid=currentResid;
      }
    return highResid;
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
    
    //if(_level==_options.maxLevels-1 && _level!=0)
    //doSweeps(10);
    //else
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

	const int cellCount=coarserCells.getSelfCount();
	const int cellTotCount=coarserCells.getCount();

	const int klen=finerKspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=finerKspace.getkvol(k);
	    Tkvol& ckvol=coarserKspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Tmode& cmode=ckvol.getmode(m);
		Field& varField=mode.getfield();
		Field& cvarField=cmode.getfield();
		Field& resField=mode.getresid();
		Field& cinjField=cmode.getInjected();
		Field& cFASField=cmode.getFASfield();
		TArray& coarserVar=dynamic_cast<TArray&>(cvarField[coarserCells]);
		TArray& coarserInj=dynamic_cast<TArray&>(cinjField[coarserCells]);
		TArray& coarserFAS=dynamic_cast<TArray&>(cFASField[coarserCells]);
		TArray& finerVar=dynamic_cast<TArray&>(varField[finerCells]);
		TArray& finerRes=dynamic_cast<TArray&>(resField[finerCells]);
		
		for(int c=0;c<cellCount;c++)
		  {
		    const int fineCount=CoarserToFiner.getCount(c);
		    coarserVar[c]=0.;
		    coarserFAS[c]=0.;
		    
		    for(int fc=0;fc<fineCount;fc++)
		      {
			const int cell=CoarserToFiner(c,fc);
			coarserVar[c]+=finerVar[cell]*finerVol[cell];
			coarserFAS[c]+=finerRes[cell];	
		      }
		    coarserVar[c]/=coarserVol[c];
		    coarserInj[c]=coarserVar[c];
		  }

		for(int c=cellCount;c<cellTotCount;c++)
		  {
		    coarserVar[c]=0.;
		    coarserFAS[c]=0.;
		    
		    const int cell=CoarserToFiner(c,0);
		    coarserVar[c]+=finerVar[cell];
		    coarserFAS[c]+=finerRes[cell];	
		    coarserInj[c]=coarserVar[c];
		  }
	       
	      }
	  }
	
	TArray& coarserVar=dynamic_cast<TArray&>(coarserMacro.temperature[coarserCells]);
	TArray& coarserInj=dynamic_cast<TArray&>(coarserMacro.TlInjected[coarserCells]);
	TArray& coarserFAS=dynamic_cast<TArray&>(coarserMacro.TlFASCorrection[coarserCells]);
	TArray& finerVar=dynamic_cast<TArray&>(_macro.temperature[finerCells]);
	TArray& finerRes=dynamic_cast<TArray&>(_macro.TlResidual[finerCells]);
	    
	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
	    coarserVar[c]=0.;
	    coarserFAS[c]=0.;
		
	    for(int fc=0;fc<fineCount;fc++)
	      {
		const int cell=CoarserToFiner(c,fc);
		coarserVar[c]+=finerVar[cell]*finerVol[cell];
		coarserFAS[c]+=finerRes[cell];
	      }
	    coarserVar[c]/=coarserVol[c];
	    coarserInj[c]=coarserVar[c];
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

	const int klen=kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Field& resField=mode.getresid();
		Field& fasField=mode.getFASfield();
		TArray& resArray=dynamic_cast<TArray&>(resField[cells]);
		TArray& fasArray=dynamic_cast<TArray&>(fasField[cells]);
		fasArray-=resArray;
	      }
	  }

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

	const int klen=kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Field& e0Field=mode.gete0field();
		TArray& e0Array=dynamic_cast<TArray&>(e0Field[cells]);

		for(int c=0;c<cellCount;c++)
		  e0Array[c]=mode.calce0(Tl[c]);

	      }
	  }
      }

    /*
    COMETInterface<T> ComInt(_meshes,_kspaces,_MeshKspaceMap,_macro,_geomFields);

    const int listSize=_IClist.size();
    for(int ic=0;ic<listSize;ic++)
      {
	COMETIC<T>* icPtr=_IClist[ic];
	const int mid0=icPtr->MeshID0;
	const int mid1=icPtr->MeshID1;
	if(icPtr->InterfaceModel=="DMM")
	  ComInt.makeDMMcoeffs(*icPtr);
	ComInt.updateOtherGhost(*icPtr,mid0);
	ComInt.updateOtherGhost(*icPtr,mid1);
      }
    */
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
	BCcellArray& BCArray=*(_BCells[n]);
	
	const int cellCount=coarserCells.getSelfCount();
	const int cellTotCount=coarserCells.getCount();
	
	const int klen=finerKspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=finerKspace.getkvol(k);
	    Tkvol& ckvol=coarserKspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Tmode& cmode=ckvol.getmode(m);
		Field& varField=mode.getfield();
		Field& cvarField=cmode.getfield();
		Field& cinjField=cmode.getInjected();
		TArray& coarserArray=dynamic_cast<TArray&>(cvarField[coarserCells]);
		TArray& finerArray=dynamic_cast<TArray&>(varField[finerCells]);
		TArray& injArray=dynamic_cast<TArray&>(cinjField[coarserCells]);
		
		for(int c=0;c<cellCount;c++)
		  {
		    const int fineCount=CoarserToFiner.getCount(c);
		    const T correction=coarserArray[c]-injArray[c];
		    
		    for(int fc=0;fc<fineCount;fc++)
		      finerArray[CoarserToFiner(c,fc)]+=correction;
		  }

		for(int c=cellCount;c<cellTotCount;c++)
		  {
		    if(BCArray[c]==4)  //correct interface cells only
		      {
			const int fineCount=CoarserToFiner.getCount(c);
			const T correction=coarserArray[c]-injArray[c];
			
			for(int fc=0;fc<fineCount;fc++)
			  finerArray[CoarserToFiner(c,fc)]+=correction;
		      }
		  }

	      }
	  }

	TArray& coarserArray=dynamic_cast<TArray&>(coarserMacro.temperature[coarserCells]);
	TArray& injArray=dynamic_cast<TArray&>(coarserMacro.TlInjected[coarserCells]);
	TArray& finerArray=dynamic_cast<TArray&>(_macro.temperature[finerCells]);
	
	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
	    const T correction=coarserArray[c]-injArray[c];
	    
	    for(int fc=0;fc<fineCount;fc++)
	      finerArray[CoarserToFiner(c,fc)]+=correction;
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
    
    applyTemperatureBoundaries();
    calcModeTemps();
  }

  T HeatFluxIntegral(const Mesh& mesh, const int faceGroupId)
  {
    T r(0.);
    bool found = false;
    const int n=mesh.getID();
    Tkspace& kspace=*_kspaces[_MeshKspaceMap[n]];
    const T DK3=kspace.getDK3();
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
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
	    
	    for(int k=0;k<kspace.getlength();k++)
	      {
		Tkvol& kv=kspace.getkvol(k);
		int modenum=kv.getmodenum();
		for(int m=0;m<modenum;m++)
		  {
		    VectorT3 vg=kv.getmode(m).getv();
		    T dk3=kv.getdk3();
		    Field& efield=kv.getmode(m).getfield();
		    const TArray& eval=dynamic_cast<const TArray&>(efield[cells]);
		    for(int f=0; f<nFaces; f++)
		      {
			const VectorT3 An=faceArea[f];
			const int c1=faceCells(f,1);
			const T vgdotAn=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
			r += eval[c1]*vgdotAn*dk3/DK3;
		      }
		  }
	      }
	    found=true;
	  }
      }
    if (!found)
      throw CException("getHeatFluxIntegral: invalid faceGroupID");
    return r*DK3;
  }

  T getWallArea(const Mesh& mesh, const int faceGroupId)
  {
    T r(0.);
    bool found = false;
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
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
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
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
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh=*_meshes[n];
	const StorageSite& cells=mesh.getCells();
	const int numcells=cells.getCount();
	const int selfcells=cells.getSelfCount();
	TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[cells]);
	const VectorT3Array& coords=dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
	T xmax(0.);
	T xmin(0.);

	for(int c=0;c<numcells;c++)
	  {
	    T x=coords[c][0];
	    if(x>xmax)
	      xmax=x;
	    if(x<xmin)
	      xmin=x;
	  }
	
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
  
  void setBCMap(COMETBCMap* bcMap) {_bcMap=*bcMap;}
  void setCoarserLevel(TCOMET* cl) {_coarserLevel=cl;}
  void setFinerLevel(TCOMET* fl) {_finerLevel=fl;}
  int getLevel() {return _level;}
  const MeshList& getMeshList() {return _meshes;}
  GeomFields& getGeomFields() {return _geomFields;}
  TkspList& getKspaces() {return _kspaces;}
  Tkspace& getKspace(const int i) {return *_kspaces[i];}
  PhononMacro& getMacro() {return _macro;}
  T getResidual() {return _residual;}

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

};

#endif
