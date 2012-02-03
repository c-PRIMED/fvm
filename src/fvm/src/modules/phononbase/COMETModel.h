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


template<class T>
class COMETModel : public Model
{

 public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
  typedef shared_ptr<VectorT3Array> VT3Ptr;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef Array<bool> BArray;
  typedef Array<int> IntArray;
  typedef Array<T> TArray;
  typedef shared_ptr<TArray> TArrptr;
  typedef map<int,COMETBC<T>*> COMETBCMap;
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
  
 COMETModel(const MeshList& meshes, const int level, GeomFields& geomFields,
	    Tkspace& kspace, PhononMacro& macro):
  Model(meshes),
    _level(level),
    _geomFields(geomFields),
    _kspace(kspace),
    _macro(macro),
    _residual(0.0)
    {
      const int numMeshes = _meshes.size();
      
      for (int n=0; n<numMeshes; n++)
	{
	  const Mesh& mesh = *_meshes[n];
	  const StorageSite& faces=mesh.getFaces();
	  const StorageSite& cells=mesh.getCells();
	  const int faceCount=faces.getCount();
	  const int cellCount=cells.getSelfCount();
	  
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
		  if (_bcMap.find(fg.id) == _bcMap.end())
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

  void init()
  {
    const int numMeshes=_meshes.size();
    const T Tinit=_options["initialTemperature"];

     for (int n=0;n<numMeshes;n++)
       {
	 const Mesh& mesh=*_meshes[n];
	 const int numK=_kspace.getlength();
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
	 
	 T e0sum=0.;
	 
	 for (int k=0;k<numK;k++)
	   {
	     Tkvol& kv=_kspace.getkvol(k);
	     const int numM=kv.getmodenum();
	     const T dk3=kv.getdk3();

	     for (int m=0;m<numM;m++)
	       {
		 Tmode& mode=kv.getmode(m);
		 T tau=mode.gettau();
		 Field& efield=mode.getfield();
		 Field& e0field=mode.gete0field();
		 Field& resfield=mode.getresid();
		 TArrptr e0var=shared_ptr<TArray>(new TArray(numcells));
		 TArrptr evar=shared_ptr<TArray>(new TArray(numcells));
		 TArrptr resid=shared_ptr<TArray>(new TArray(numcells));
		 const T einit=mode.calce0(Tinit);
		 e0sum+=einit*dk3/tau;
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
		 
		 BCfaceArray& BCfArray=*(_BFaces[n]);
		 BCcellArray& BCArray=*(_BCells[n]);
		 foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		   {
		     const FaceGroup& fg = *fgPtr;
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
			     Refl_pair refls;
			     Refl_pair reflsFrom;
			     _kspace.findSpecs(dk3,vmag,m,so,refls);
			     rmap[fg.id]=refls;
			     const int k1=refls.first.second;
			     Tmode& mode2=_kspace.getkvol(k1).getmode(m);
			     Refl_Map& rmap2=mode2.getreflmap();
			     reflsFrom.first.second=-1;
			     reflsFrom.second.second=k;
			     rmap2[fg.id]=reflsFrom;
			   }
		       }
		     else
		       {
			 const StorageSite& faces = fg.site;
			 const int faceCount=faces.getCount();
			 const int offSet=faces.getOffset();
			 Refl_pair refls;
			 refls.first.second=-1;
			 refls.second.second=-1;
			 rmap[fg.id]=refls;

			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=1;
		       }
		   }
	       }    
	   }
	 
	 shared_ptr<TArray> TlResidCell(new TArray(numcells));
	 *TlResidCell=0.;
	 _macro.TlResidual.addArray(cells,TlResidCell);
       }
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
	 const Mesh& mesh=*_meshes[n];
	 const int numK=_kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 shared_ptr<TArray> deltaTcell(new TArray(numcells));
	 *deltaTcell=0.; 
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 _macro.deltaT.addArray(cells,deltaTcell);
	 
	 T e0sum=0.;
	 
	 for (int k=0;k<numK;k++)
	   {
	     Tkvol& kv=_kspace.getkvol(k);
	     const int numM=kv.getmodenum();
	     const T dk3=kv.getdk3();

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
		 e0sum+=einit*dk3;
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

		 BCfaceArray& BCfArray=*(_BFaces[n]);
		 BCcellArray& BCArray=*(_BCells[n]);
		 foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
		   {
		     const FaceGroup& fg = *fgPtr;
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
		     else
		       {
			 const StorageSite& faces = fg.site;
			 const int faceCount=faces.getCount();
			 const int offSet=faces.getOffset();
			 
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=1;
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
     applyTemperatureBoundaries();
  }

  void applyTemperatureBoundaries()
  {
    const int numMeshes=_meshes.size();

    for(int n=0;n<numMeshes;n++)
      {
	const Mesh& mesh=*_meshes[n];
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg = *fgPtr;
	    const StorageSite& faces = fg.site;
	    const COMETBC<T>& bc = *_bcMap[fg.id];
	    
	    COMETBoundary<T> cbc(faces, mesh,_geomFields,_kspace,_options,fg.id);

	    if(bc.bcType=="temperature")
	      {	      
		FloatValEvaluator<T>
		  bTemperature(bc.getVal("specifiedTemperature"),faces);
	  
		cbc.applyTemperatureWall(bTemperature);
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
	    Tkspace* newKspacePtr=new Tkspace();
	    PhononMacro* newMacroPtr=new PhononMacro("coarse");

	    newKspacePtr->CopyKspace(finerModel->getKspace());

	    int newCount= MakeCoarseMesh(finerModel->getMeshList(),
					 finerModel->getGeomFields(),
					 *newMeshesPtr);
	    	    
	    TCOMET* newModelPtr=new COMETModel(*newMeshesPtr,thisLevel,
					       finerModel->getGeomFields(),
					       *newKspacePtr,*newMacroPtr);
	    
	    newModelPtr->setFinerLevel(finerModel);
	    finerModel->setCoarserLevel(newModelPtr);
	    newModelPtr->getOptions()=finerModel->getOptions();
	    newModelPtr->getBCs()=finerModel->getBCs();

	    newModelPtr->initCoarse();
	    newModelPtr->MakeCoarseModel(newModelPtr);
	  }
      }
    else if(_options.AgglomerationMethod=="AMG")
      throw CException("Have not implemented AMG agglomeration method.");
    else
      throw CException("Unknown agglomeration method.");
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
	foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
	  {
	    const FaceGroup& fg=*fgPtr;
	    const int size=fg.site.getCount();
	    newMeshPtr->createBoundaryFaceGroup(size,inOffset,fg.id,fg.groupType);
	    inOffset+=size;
	  }

	//now make the geom fields
	const int outCellsCount=outCells.getSelfCount();
	TArrptr outCellVolumePtr=TArrptr(new TArray(outCellsCount));
	TArray& outCV=*outCellVolumePtr;
	outCV=0.;

	Field& VolumeField=inGeomFields.volume;
	const TArray& inCV=dynamic_cast<const TArray&>(VolumeField[inCells]);

	for(int c=0;c<outCellsCount;c++)
	  {
	    const int fineCount=CoarseToFineCells->getCount(c);
	    for(int i=0;i<fineCount;i++)
	      {
		int fc=(*CoarseToFineCells)(c,i);
		outCV[c]+=inCV[fc];
	      }
	  }
	
	VolumeField.addArray(outCells,outCellVolumePtr);

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
      smooth();
  }

  void smooth()
  {
    const T relFactor=_options.relFactor;
    const int numMeshes=_meshes.size();
    for(int msh=0;msh<numMeshes;msh++)
      {
	const Mesh& mesh=*_meshes[msh];
	const BCcellArray& BCArray=*(_BCells[msh]);
	const BCfaceArray& BCfArray=*(_BFaces[msh]);
	COMETDiscretizer<T> CDisc(mesh,_geomFields,_macro,
				  _kspace,_bcMap,BCArray,BCfArray,_options);

	CDisc.setfgFinder();
	CDisc.COMETSolve(1,_level); //forward
	CDisc.COMETSolve(-1,_level); //reverse
      }
  }

  T updateResid(const bool addFAS)
  {
    const int numMeshes=_meshes.size();
    T lowResid=-1.;
    T currentResid;
    for(int msh=0;msh<numMeshes;msh++)
      {
	const Mesh& mesh=*_meshes[msh];
	const BCcellArray& BCArray=*(_BCells[msh]);
	const BCfaceArray& BCfArray=*(_BFaces[msh]);
	COMETDiscretizer<T> CDisc(mesh,_geomFields,_macro,
				  _kspace,_bcMap,BCArray,BCfArray,_options);
	
	CDisc.setfgFinder();
	CDisc.findResid(addFAS);
	currentResid=CDisc.getAveResid();

	if(lowResid<0)
	  lowResid=currentResid;
	else
	  if(currentResid<lowResid)
	    lowResid=currentResid;
      }
    return lowResid;
  }

  void cycle()
  {
    doSweeps(_options.preSweeps);
    
    if(_level+1<_options.maxLevels)
      {
	if(_level==0)
	  updateResid(false);
	else
	  updateResid(true);

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
	const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
	Tkspace& coarserKspace=_coarserLevel->getKspace();
	PhononMacro& coarserMacro=_coarserLevel->getMacro();
	const StorageSite& finerCells=finerMesh.getCells();
	const StorageSite& coarserCells=coarserMesh.getCells();
	const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);
	const TArray& coarserVol=dynamic_cast<TArray&>(_geomFields.volume[coarserCells]);
	const TArray& finerVol=dynamic_cast<TArray&>(_geomFields.volume[finerCells]);

	const int cellCount=coarserCells.getSelfCount();

	const int klen=_kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
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
	const StorageSite& cells=mesh.getCells();

	const int klen=_kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
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
	const StorageSite& cells=mesh.getCells();
	const int cellCount=cells.getSelfCount();
	TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[cells]);

	const int klen=_kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
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
  }

  void correctSolution()
  {
    const int numMeshes = _meshes.size();
    
    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& finerMesh=*_meshes[n];
	const Mesh& coarserMesh=*(_coarserLevel->getMeshList())[n];
	Tkspace& coarserKspace=_coarserLevel->getKspace();
	PhononMacro& coarserMacro=_coarserLevel->getMacro();
	const StorageSite& finerCells=finerMesh.getCells();
	const StorageSite& coarserCells=coarserMesh.getCells();
	const CRConnectivity& CoarserToFiner=coarserMesh.getConnectivity(coarserCells,finerCells);
	
	const int cellCount=coarserCells.getSelfCount();
	
	const int klen=_kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
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
	const StorageSite& cells = mesh.getCells();
	const int numcells = cells.getCount();
	const T tauTot=_kspace.calcTauTot();
	TArray& TL=dynamic_cast<TArray&>(_macro.temperature[cells]);
	TArray& e0Array=dynamic_cast<TArray&>(_macro.e0[cells]);
	
	for(int c=0;c<numcells;c++)
	  _kspace.NewtonSolve(TL[c],e0Array[c]*tauTot);
      }
  }
  
  void advance(const int iters)
  {
    applyTemperatureBoundaries();
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
    //cout<<endl<<"Total Iterations:"<<niters<<" Residual:"<<_residual<<endl;
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
		    const TArray& eval=dynamic_cast<const TArray&>(efield[cells]);
		    for(int f=0; f<nFaces; f++)
		      {
			const VectorT3 An=faceArea[f];
			const int c1=faceCells(f,1);
			const T vgdotAn=An[0]*vg[0]+An[1]*vg[1]+An[2]*vg[2];
			r += eval[c1]*vgdotAn*dk3;
		      }
		  }
	      }
	    found=true;
	  }
      }
    if (!found)
      throw CException("getHeatFluxIntegral: invalid faceGroupID");
    return r;
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
    if (!found)
      throw CException("getwallArea: invalid faceGroupID");
    return An;
  }
  
  ArrayBase* getValueArray(const Mesh& mesh, const int cell)
    {
      //only returns the e" values, not the lattice temperature
      const int allModes=_kspace.gettotmodes();
      TArray* vals=new TArray(allModes);
      const StorageSite& cells=mesh.getCells();
      const int len=_kspace.getlength();
      int count=0;
      for(int k=0;k<len;k++)
	{
	  Tkvol& kvol=_kspace.getkvol(k);
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
  
  void setBCMap(COMETBCMap* bcMap) {_bcMap=*bcMap;}
  void setCoarserLevel(TCOMET* cl) {_coarserLevel=cl;}
  void setFinerLevel(TCOMET* fl) {_finerLevel=fl;}
  int getLevel() {return _level;}
  const MeshList& getMeshList() {return _meshes;}
  GeomFields& getGeomFields() {return _geomFields;}
  Tkspace& getKspace() {return _kspace;}
  PhononMacro& getMacro() {return _macro;}
  T getResidual() {return _residual;}

 private:

  const int _level;    //0 being the finest level
  GeomFields& _geomFields;
  Tkspace& _kspace;
  PhononMacro& _macro;
  TCOMET* _finestLevel;
  TCOMET* _coarserLevel;
  TCOMET* _finerLevel;
  COMETBCMap _bcMap;
  TCModOpts _options;
  BCfaceList _BFaces;
  BCcellList _BCells;
  T _residual;

};

#endif
