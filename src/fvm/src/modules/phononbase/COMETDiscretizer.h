#ifndef _COMETDISCRETIZER_H_
#define _COMETDISCRETIZER_H_

#include "Mesh.h"
#include "Kspace.h"
#include "kvol.h"
#include "pmode.h"
#include "ArrowHeadMatrix.h"
#include "Array.h"
#include "Vector.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "GeomFields.h"
#include "NumType.h"
#include "COMETBC.h"
#include <math.h>
#include <map>
#include "MatrixJML.h"
#include "SquareMatrix.h"
#include "PhononMacro.h"
#include "SquareTensor.h"

template<class T>
class COMETDiscretizer
{

 public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef Vector<T_Scalar,3> VectorT3;  
  typedef Array<VectorT3> VectorT3Array;
  typedef Array<T_Scalar> TArray;
  typedef MatrixJML<T> TMatrix;
  typedef ArrowHeadMatrix<T> TArrow;
  typedef SquareMatrix<T> TSquare;
  typedef map<int,COMETBC<T>*> COMETBCMap;
  typedef COMETModelOptions<T> COpts;
  typedef Array<int> IntArray;
  typedef Array<bool> BoolArray;
  typedef Vector<int,2> VecInt2;
  typedef map<int,VecInt2> FaceToFg;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef SquareTensor<T,3> T3Tensor;

 COMETDiscretizer(const Mesh& mesh, const GeomFields& geomfields, 
		  PhononMacro& macro, Tkspace& kspace, COMETBCMap& bcMap,
		  const IntArray& BCArray, const IntArray& BCfArray, COpts& options):
  _mesh(mesh),
    _geomFields(geomfields),
    _cells(mesh.getCells()),
    _faces(mesh.getFaces()),
    _cellFaces(mesh.getCellFaces()),
    _faceCells(mesh.getAllFaceCells()),
    _areaMagField(_geomFields.areaMag),
    _faceAreaMag((dynamic_cast<const TArray&>(_areaMagField[_faces]))),
    _areaField(_geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _cellVolume(dynamic_cast<const TArray&>(_geomFields.volume[_cells])),
    _macro(macro),
    _kspace(kspace),
    _bcMap(bcMap),
    _BCArray(BCArray),
    _BCfArray(BCfArray),
    _aveResid(-1.),
    _residChange(-1.),
    _fgFinder(),
    _options(options)
    {}

  void COMETSolve(const int dir,const int level)
  {
    const int cellcount=_cells.getSelfCount();
    int start;

    if(dir==1)
      start=0;
    if(dir==-1)
      start=cellcount-1;
    
    for(int c=start;((c<cellcount)&&(c>-1));c+=dir)
      {	
	const int totalmodes=_kspace.gettotmodes();
	TArray Bvec(totalmodes+1);
	TArray Resid(totalmodes+1);
	
	if(_BCArray[c]==0)  //not a reflecting boundary
	  {
	    TArrow AMat(totalmodes+1);
	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    COMETConvection(c,AMat,Bvec);
	    COMETCollision(c,&AMat,Bvec);
	    COMETEquilibrium(c,&AMat,Bvec);

	    if(_options.withNormal)
	      COMETShifted(c,&AMat,Bvec);
	    
	    if(level>0)
	      addFAS(c,Bvec);
	    
	    Resid=Bvec;
	    AMat.Solve(Bvec);
	    
	    Distribute(c,Bvec,Resid);
	  }
	else if(_BCArray[c]==1) //reflecting boundary
	  {
	    for(int j=0;j<3;j++)
	      {
		TSquare AMat(totalmodes+1);
		Bvec.zero();
		Resid.zero();
		AMat.zero();
		
		COMETConvection(c,AMat,Bvec);	
		COMETCollision(c,&AMat,Bvec);
		COMETEquilibrium(c,&AMat,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
		
		Resid=Bvec;
		AMat.Solve(Bvec);
		
		Distribute(c,Bvec,Resid);
		updatee0();
	      }
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");

      }
    updatee0();
    if(_options.withNormal)
      updateeShifted();
  }

  void COMETConvection(const int cell, TArrow& Amat, TArray& BVec)
  {
    /* This is the COMET discretization for cells that do not
       do not have a face which is a reflecting boundary.  When
       there is a face with a reflecting boundary, we can no longer
       use an arrowhead matrix as the structure becomes unknown a priori.
     */
    const int neibcount=_cellFaces.getCount(cell);
    
    for(int j=0;j<neibcount;j++)
      {
	const int f=_cellFaces(cell,j);
	int cell2=_faceCells(f,1);
	const int klen=_kspace.getlength();
	VectorT3 Af=_faceArea[f];

	T flux;
	int count=1;

	if(cell2==cell)
	  {
	    Af=Af*(-1.);
	    cell2=_faceCells(f,0);
	  }
	
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		VectorT3 vg=mode.getv();
		Field& efield=mode.getfield();
		TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
		flux=vg[0]*Af[0]+vg[1]*Af[1]+vg[2]*Af[2];
		if(flux>T_Scalar(0))
		  {
		    Amat.getElement(count,count)-=flux;
		    BVec[count-1]-=flux*eArray[cell];
		  }
		else
		  BVec[count-1]-=flux*eArray[cell2];
		
		count+=1;
	      }
	  }	   	
      }   
  }

  void COMETConvection(const int cell, TSquare& Amat, TArray& BVec)
  {
    const int neibcount=_cellFaces.getCount(cell);

    for(int j=0;j<neibcount;j++)
      {
	const int f=_cellFaces(cell,j);
	int cell2=_faceCells(f,1);
	const int klen=_kspace.getlength();
	VectorT3 Af=_faceArea[f];
	
	T flux;
	int count=1;

	if(cell2==cell)
	  {
	    Af=Af*(-1.);
	    cell2=_faceCells(f,0);
	  }
	
	if(_BCfArray[f]==2) //If the face in question is a reflecting face
	  {
	    int Fgid=findFgId(f);
	    T refl=(*(_bcMap[Fgid]))["specifiedReflection"];
	    T oneMinusRefl=1-refl;

	    //first sweep - have to make sumVdotA
	    T sumVdotA=0.;
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		T dk3=kvol.getdk3();
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    Field& efield=mode.getfield();
		    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
		    flux=vg[0]*Af[0]+vg[1]*Af[1]+vg[2]*Af[2];
		    if(flux>T_Scalar(0))
		      {
			Amat(count,count)-=flux;
			BVec[count-1]-=flux*eArray[cell];
			sumVdotA+=flux*dk3;
		      }
		    else
		      {//have to move through all other modes
			int ccount=1;
			for(int kk=0;kk<klen;kk++)
			  {
			    Tkvol& kkvol=_kspace.getkvol(kk);
			    T ddk3=kkvol.getdk3();
			    const int numMODES=kkvol.getmodenum();
			    for(int mm=0;mm<numMODES;mm++)
			      {
				Refl_pair& rpairs=kkvol.getmode(mm).getReflpair(Fgid);
				int k1=rpairs.first.second;
				if(kk!=k || m!=mm)
				  {
				    VectorT3 vvg=kkvol.getmode(mm).getv();
				    T VdotA=vvg[0]*Af[0]+vvg[1]*Af[1]+vvg[2]*Af[2];    
				    if(VdotA>T_Scalar(0))
				      {
					if(k==k1 && m==mm)
					  Amat(count,ccount)-=VdotA*refl*ddk3/dk3;
				      }
				  }
				ccount++;
			      }
			  }
		      }
		    count++;
		  }
	      }

	    //Second sweep
	    count=1;
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		T dk3=kvol.getdk3();
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    VectorT3 vg=kvol.getmode(m).getv();
		    flux=vg[0]*Af[0]+vg[1]*Af[1]+vg[2]*Af[2];
		    if(flux<T_Scalar(0))
		      {
	       		int ccount=1;
			for(int kk=0;kk<klen;kk++)
			  {
			    Tkvol& kkvol=_kspace.getkvol(kk);
			    const int numMODES=kkvol.getmodenum();
			    T ddk3=kkvol.getdk3();
			    for(int mm=0;mm<numMODES;mm++)
			      {
				if(kk!=k || m!=mm)
				  {
				    VectorT3 vvg=kkvol.getmode(mm).getv();
				    T VdotA=vvg[0]*Af[0]+vvg[1]*Af[1]+vvg[2]*Af[2];
				    if(VdotA>T_Scalar(0))
				      {
					Refl_pair& rpairs=kkvol.getmode(mm).getReflpair(Fgid);
					int k1=rpairs.first.second;
					Field& eefield=kvol.getmode(mm).getfield();
					TArray& eeArray=dynamic_cast<TArray&>(eefield[_cells]);
					Amat(count,ccount)-=flux*VdotA*ddk3*oneMinusRefl/sumVdotA;
					BVec[count-1]-=flux*VdotA*ddk3
					  *eeArray[cell]/sumVdotA*oneMinusRefl;
					if(k==k1 && m==mm)
					  BVec[count-1]-=eeArray[cell]*ddk3/dk3*refl*VdotA;
				      }
				  }
				ccount++;
			      }
			  }
		      }
		    count++;
		  }
	      }


	  }
	else  //if the face in question is not reflecting
	  {
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    Field& efield=mode.getfield();
		    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
		    flux=vg[0]*Af[0]+vg[1]*Af[1]+vg[2]*Af[2];
		    if(flux>T_Scalar(0))
		      {
			Amat(count,count)-=flux;
			BVec[count-1]-=flux*eArray[cell];
		      }
		    else
		      BVec[count-1]-=flux*eArray[cell2];
		    
		    count+=1;
		  }
	      }
	  }
	
	
      }
    
  }
  
  void COMETCollision(const int cell, TMatrix* Amat, TArray& BVec)
  {
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    TArray& Tlold=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    int count=1;
    T coeff;
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    T tau=mode.gettau();
	    T de0dT=mode.calcde0dT(Tlold[cell]);
	    Field& efield=mode.getfield();
	    Field& e0field=mode.gete0field();
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
	    TArray& e0Array=dynamic_cast<TArray&>(e0field[_cells]);
	    coeff=_cellVolume[cell]/tau;
	    Amat->getElement(count,order)+=coeff*de0dT;
	    Amat->getElement(count,count)-=coeff;
	    BVec[count-1]-=coeff*eArray[cell];
	    //BVec[count-1]+=coeff*Tlold[cell]*de0dT;
	    BVec[count-1]+=coeff*e0Array[cell];
	    count+=1;
	  }
      }
    
  }

  void COMETEquilibrium(const int cell, TMatrix* Amat, TArray& BVec)
  {
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    TArray& Tlold=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    const T tauTot=_kspace.getde0taudT(Tlold[cell]);
    int count=1;
    T coeff;
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    T tau=mode.gettau();
	    coeff=dk3/tau/tauTot;
	    Field& efield=mode.getfield();
	    Field& e0field=mode.gete0field();
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
	    TArray& e0Array=dynamic_cast<TArray&>(e0field[_cells]);
	    Amat->getElement(order,count)+=coeff;
	    BVec[totalmodes]+=coeff*eArray[cell];
	    BVec[totalmodes]-=coeff*e0Array[cell];
	    count+=1;
	  }
      }
    Amat->getElement(order,order)=-1.;
    //BVec[totalmodes]+=Tlold[cell];
  }

  void COMETShifted(const int cell, TMatrix* Amat, TArray& BVec)
  { //adds to collision and equilibrium
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    TArray& Tlold=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    const T tauTot=_kspace.getde0taudT(Tlold[cell]);
    int count=1;
    T coeff;
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    T tau=mode.gettauN();
	    Field& efield=mode.getfield();
	    Field& eShiftedfield=mode.geteShifted();
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
	    TArray& eShifted=dynamic_cast<TArray&>(eShiftedfield[_cells]);
	    coeff=_cellVolume[cell]/tau;
	    Amat->getElement(count,count)-=coeff;
	    Amat->getElement(order,count)+=coeff/tauTot;
	    BVec[count-1]-=coeff*eArray[cell];
	    BVec[count-1]+=coeff*eShifted[cell];
	    BVec[totalmodes]+=coeff*eArray[cell]/tauTot;
	    BVec[totalmodes]-=coeff*eShifted[cell]/tauTot;
	    count+=1;
	  }
      }
    
  }

  void Distribute(const int cell, TArray& BVec, TArray& Rvec)
  {
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    int count=1;

    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Field& efield=mode.getfield();
	    Field& resfield=mode.getresid();
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
	    TArray& resArray=dynamic_cast<TArray&>(resfield[_cells]);
	    eArray[cell]-=BVec[count-1];
	    resArray[cell]=-Rvec[count-1];
	    count+=1;
	  }
      }
    
    TArray& TlArray=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    TlArray[cell]-=BVec[totalmodes];
    TArray& deltaTArray=dynamic_cast<TArray&>(_macro.deltaT[_cells]);
    deltaTArray[cell]=BVec[totalmodes];
    TArray& e0ResArray=dynamic_cast<TArray&>(_macro.e0Residual[_cells]);
    e0ResArray[cell]=-Rvec[totalmodes];
  }

  void findResid(const bool plusFAS)
  {
    const int cellcount=_cells.getSelfCount();
    const int totalmodes=_kspace.gettotmodes();
    TArray ResidSum(totalmodes+1);
    TArray Bsum(totalmodes+1);
    T ResidScalar=0.;
    T traceSum=0.;
    TArray Bvec(totalmodes+1);
    TArray Resid(totalmodes+1);
    ResidSum.zero();
    Bsum.zero();
    
    for(int c=0;c<cellcount;c++)
      {	
	if(_BCArray[c]==0)  //not a reflecting boundary
	  {
	    TArrow AMat(totalmodes+1);
	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    COMETConvection(c,AMat,Bvec);	
	    COMETCollision(c,&AMat,Bvec);
	    COMETEquilibrium(c,&AMat,Bvec);
	    
	    if(plusFAS)
	      addFAS(c,Bvec);

	    traceSum+=AMat.getTraceAbs();
	    Resid+=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);

	    makeValueArray(c,Bvec);
	    ArrayAbs(Bvec);
	    ArrayAbs(Resid);
	    Bsum+=Bvec;
	    ResidSum+=Resid;	    
	  }
	else if(_BCArray[c]==1) //reflecting boundary
	  {
	    TSquare AMat(totalmodes+1);
	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    COMETConvection(c,AMat,Bvec);	
	    COMETCollision(c,&AMat,Bvec);
	    COMETEquilibrium(c,&AMat,Bvec);
	    
	    if(plusFAS)
	      addFAS(c,Bvec);

	    traceSum+=AMat.getTraceAbs();
	    Resid+=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);

	    makeValueArray(c,Bvec);
	    ArrayAbs(Resid);
	    ArrayAbs(Bvec);
	    Bsum+=Bvec;
	    ResidSum+=Resid;
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");
      }

    for(int o=0;o<totalmodes+1;o++)
      {
	ResidScalar+=ResidSum[o];
      }

    ResidScalar/=traceSum;

    if(_aveResid==-1)
      {_aveResid=ResidScalar;}
    else
      {
	_residChange=fabs(_aveResid-ResidScalar)/_aveResid;
	_aveResid=ResidScalar;
      }
  }
  
  TArray gatherResid(const int c)
  {
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    int count=1;

    TArray rArray(order);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Field& resfield=mode.getresid();
	    TArray& resArray=dynamic_cast<TArray&>(resfield[_cells]);
	    rArray[count]=resArray[c];
	    count+=1;
	  }
      }
    rArray[order]=_macro.e0Residual[c];
  }
  
  T getResidChange() {return _residChange;}
  T getAveResid() {return _aveResid;}

  void setfgFinder()
  {
    foreach(const FaceGroupPtr fgPtr, _mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg=*fgPtr;
	const int off=fg.site.getOffset();
	const int cnt=fg.site.getCount();
	const int id=fg.id;
	VecInt2 BegEnd;
	
	BegEnd[0]=off;
	BegEnd[1]=off+cnt;
	_fgFinder[id]=BegEnd;
      }
  }

  int findFgId(const int faceIndex)
  {
    FaceToFg::iterator id;
    for(id=_fgFinder.begin();id!=_fgFinder.end();id++)
      {
	if(id->second[0]<=faceIndex && id->second[1]>faceIndex)
	  return id->first;
      }
    throw CException("Didn't find matching FaceGroup!");
    return -1;
  }

  void ArrayAbs(TArray& o)
  {
    int length=o.getLength();
    for(int i=0;i<length;i++)
      o[i]=fabs(o[i]);
  }

  void makeValueArray(const int c, TArray& o)
  {
    int klen=_kspace.getlength();
    int count=0;
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Field& efield=mode.getfield();
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);
	    o[count]=eArray[c];
	    count+=1;
	  }
      }
  }

  void addFAS(const int c, TArray& bVec)
  {
    int klen=_kspace.getlength();
    int count=0;
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Field& fasField=mode.getFASfield();
	    TArray& fasArray=dynamic_cast<TArray&>(fasField[_cells]);
	    bVec[count]-=fasArray[c];
	    count+=1;
	  }
      }
    TArray& fasArray=dynamic_cast<TArray&>(_macro.e0FASCorrection[_cells]);
    bVec[count]-=fasArray[c];
  }

  void updatee0()
  {
    TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    TArray& dT=dynamic_cast<TArray&>(_macro.deltaT[_cells]);
    const T cellCount=_cells.getSelfCount();
    int klen=_kspace.getlength();
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    Field& e0Field=mode.gete0field();
	    TArray& e0Array=dynamic_cast<TArray&>(e0Field[_cells]);
	    
	    for(int c=0;c<cellCount;c++)
	      {
		T Told=Tl[c]-dT[c];
		T de0dTold=mode.calcde0dT(Told);
		e0Array[c]-=dT[c]*de0dTold;
	      }
	  }
      }
  }

  void updateeShifted()
  {
    const int cellcount=_cells.getSelfCount();
    TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    VectorT3Array& lam=dynamic_cast<VectorT3Array&>(_macro.lam[_cells]);

    for(int c=0;c<cellcount;c++)
      {	
	T3Tensor TensorSum;
	VectorT3 VectorSum;
	VectorT3 lambda;
	T TpreFactor=0.;
	T VpreFactor=0.;

	int klen=_kspace.getlength();
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Field& eField=mode.getfield();
		TArray& eArray=dynamic_cast<TArray&>(eField[_cells]);
		TpreFactor+=mode.calcTensorPrefactor(Tl[c]);
		VpreFactor+=mode.calcVectorPrefactor(Tl[c],eArray[c]);
	      }
	  }

	TensorSum.zero();
	VectorSum.zero();

	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    T dk3=kvol.getdk3();
	    VectorT3 Kvec=kvol.getkvec();
	    T3Tensor TempTensor;
	    outerProduct(Kvec,Kvec,TempTensor);
	    TensorSum+=TempTensor*dk3*TpreFactor;
	    VectorSum+=Kvec*dk3*VpreFactor;
	  }

	lambda=inverse(TensorSum)*VectorSum;
	lam[c]=lambda;

	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    VectorT3 Kvec=kvol.getkvec();
	    T shift=Kvec[0]*lambda[0]+Kvec[1]*lambda[1]+Kvec[2]*lambda[2];
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		Field& eShiftedField=mode.geteShifted();
		TArray& eShifted=dynamic_cast<TArray&>(eShiftedField[_cells]);
		eShifted[c]=mode.calcShifted(Tl[c],shift);
	      }
	  }
	
      }
  }
  
  void outerProduct(const VectorT3& v1, const VectorT3& v2, T3Tensor& out)
  {
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	out(i,j)=v1[i]*v2[j];
  }
  
 private:
  const Mesh& _mesh;
  const GeomFields& _geomFields;
  const StorageSite& _cells;
  const StorageSite& _faces;
  const CRConnectivity& _cellFaces;
  const CRConnectivity& _faceCells;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  const TArray& _cellVolume;
  PhononMacro& _macro;
  Tkspace& _kspace;
  COMETBCMap& _bcMap;
  const IntArray& _BCArray;
  const IntArray& _BCfArray;
  T _aveResid;
  T _residChange;
  FaceToFg _fgFinder;
  COpts _options;
  
};


#endif
