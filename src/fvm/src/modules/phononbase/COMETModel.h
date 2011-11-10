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
    _macro(macro)
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
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 
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
		 Field& e0field=mode.gete0field();
		 Field& resfield=mode.getresid();
		 TArrptr e0var=shared_ptr<TArray>(new TArray(numcells));
		 TArrptr evar=shared_ptr<TArray>(new TArray(numcells));
		 TArrptr resid=shared_ptr<TArray>(new TArray(numcells));
		 const T einit=mode.calce0(Tinit);
		 e0sum+=einit*dk3;
		 *evar=einit;
		 *e0var=einit;
		 *resid=0.;
		 efield.addArray(cells,evar);
		 e0field.addArray(cells,e0var);
		 resfield.addArray(cells,resid);
		 
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

			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=2;

			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     BCArray[cell1]=1;
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
			     Refl_pair refls;
			     _kspace.findSpecs(dk3,vmag,m,so,refls);
			     rmap[fg.id]=refls;
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
	 
	 e0sum=e0sum/_kspace.getDK3();
	 shared_ptr<TArray> e0cell(new TArray(numcells));
	 *e0cell=e0sum;
	 _macro.e0.addArray(cells,e0cell);
	 shared_ptr<TArray> e0ResidCell(new TArray(numcells));
	 *e0ResidCell=0.;
	 _macro.e0Residual.addArray(cells,e0ResidCell);
       }

     MakeCoarseModel(this);
  }

  void initCoarse()
  {
    const int numMeshes=_meshes.size();
    const T Tinit=0.;

     for (int n=0;n<numMeshes;n++)
       {
	 const Mesh& mesh=*_meshes[n];
	 const int numK=_kspace.getlength();
	 const StorageSite& cells=mesh.getCells();
	 const int numcells=cells.getCount();
	 
	 shared_ptr<TArray> TLcell(new TArray(numcells));
	 *TLcell=Tinit;
	 _macro.temperature.addArray(cells,TLcell);
	 
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
		 Field& e0field=mode.gete0field();
		 Field& resfield=mode.getresid();
		 Field& FASfield=mode.getFASfield();
		 TArrptr e0var=TArrptr(new TArray(numcells));
		 TArrptr evar=TArrptr(new TArray(numcells));
		 TArrptr resid=TArrptr(new TArray(numcells));
		 TArrptr FAS=TArrptr(new TArray(numcells));
		 const T einit=0.;
		 e0sum+=einit*dk3;
		 *evar=einit;
		 *e0var=einit;
		 *resid=0.;
		 *FAS=0;
		 efield.addArray(cells,evar);
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
			 
			 for(int i=offSet;i<offSet+faceCount;i++)
			   BCfArray[i]=2;

			 for(int i=0;i<faceCount;i++)
			   {
			     int cell1=BfaceCells(i,0);
			     BCArray[cell1]=1;
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
	 
	 e0sum=e0sum/_kspace.getDK3();
	 TArrptr e0cell(new TArray(numcells));
	 *e0cell=e0sum;
	 _macro.e0.addArray(cells,e0cell);
	 TArrptr e0ResidCell(new TArray(numcells));
	 *e0ResidCell=0.;
	 _macro.e0Residual.addArray(cells,e0ResidCell);
	 TArrptr e0Inj(new TArray(numcells));
	 *e0Inj=0.;
	 _macro.e0Injected.addArray(cells,e0Inj);
	 TArrptr e0FAS(new TArray(numcells));
	 *e0FAS=0.;
	 _macro.e0FASCorrection.addArray(cells,e0FAS);
       }
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
	    
	    MakeCoarseModel(newModelPtr);
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
		      if(areaMagArray[c2]>maxArea)
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
			
			if(areaMagArray[c2]>maxArea)
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
    const int numMeshes=_meshes.size();
    for(int msh=0;msh<numMeshes;msh++)
      {
	const Mesh& mesh=*_meshes[msh];
	const BCcellArray& BCArray=*(_BCells[msh]);
	const BCfaceArray& BCfArray=*(_BFaces[msh]);
	COMETDiscretizer<T> CDisc(mesh,_geomFields,_macro,
				  _kspace,_bcMap,BCArray,BCfArray);
	
	CDisc.setfgFinder();
	applyTemperatureBoundaries();
	CDisc.COMETSolve(1,_level); //forward
	applyTemperatureBoundaries();
	CDisc.COMETSolve(-1,_level); //reverse
      }
  }

  T updateResid()
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
				  _kspace,_bcMap,BCArray,BCfArray);
	
	CDisc.setfgFinder();
	applyTemperatureBoundaries();
	CDisc.findResid();
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
	if(_options.preSweeps==0)
	  updateResid();
	injectResid();
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
		Field& cinjField=cmode.gete0field();
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

	TArray& coarserVar=dynamic_cast<TArray&>(coarserMacro.e0[coarserCells]);
	TArray& coarserInj=dynamic_cast<TArray&>(coarserMacro.e0Injected[coarserCells]);
	TArray& coarserFAS=dynamic_cast<TArray&>(coarserMacro.e0FASCorrection[coarserCells]);
	TArray& finerVar=dynamic_cast<TArray&>(_macro.e0[finerCells]);
	TArray& finerRes=dynamic_cast<TArray&>(_macro.e0Residual[finerCells]);

	for(int c=0;c<cellCount;c++)
	  {
	    const int fineCount=CoarserToFiner.getCount(c);
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
    updateResid();

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

	TArray& resArray=dynamic_cast<TArray&>(_macro.e0Residual[cells]);
	TArray& fasArray=dynamic_cast<TArray&>(_macro.e0FASCorrection[cells]);
	fasArray-=resArray;
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
		Field& cinjField=cmode.gete0field();
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

	TArray& coarserArray=dynamic_cast<TArray&>(coarserMacro.e0[coarserCells]);
	TArray& injArray=dynamic_cast<TArray&>(coarserMacro.e0Injected[coarserCells]);
	TArray& finerArray=dynamic_cast<TArray&>(_macro.e0[finerCells]);
	
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
    T residual;
    applyTemperatureBoundaries();
    residual=updateResid();
    int niters=0;
    const T absTol=_options.absTolerance;
    const int show=_options.showResidual;

    while((niters<iters) && (residual>absTol))
      {
	cycle();
	applyTemperatureBoundaries();
	niters++;
	residual=updateResid();
	if(niters%show==0)
	  cout<<"Iteration:"<<niters<<" Residual:"<<residual<<endl;
	
      }
    updateTL();
    cout<<"Total Iterations:"<<niters<<" Residual:"<<residual<<endl;
  }
  
  void setBCMap(COMETBCMap* bcMap) {_bcMap=*bcMap;}
  void setCoarserLevel(TCOMET* cl) {_coarserLevel=cl;}
  void setFinerLevel(TCOMET* fl) {_finerLevel=fl;}
  int getLevel() {return _level;}
  const MeshList& getMeshList() {return _meshes;}
  GeomFields& getGeomFields() {return _geomFields;}
  Tkspace& getKspace() {return _kspace;}
  PhononMacro& getMacro() {return _macro;}

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

};

#endif
