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
#include "PhononBC.h"
#include <math.h>
#include <map>
#include "MatrixJML.h"
#include "SquareMatrix.h"

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
  typedef ArrowHeadMatrix<T> TArrow;
  typedef SquareMatrix<T> TSquare;
  typedef map<int,PhononBC<T>*> PhononBCMap;
  typedef Array<int> IntArray;
  typedef Array<bool> BoolArray;
  typedef Vector<int,2> VecInt2;
  typedef map<int,VecInt2> FaceToFg;

 COMETDiscretizer(const Mesh& mesh, const GeomFields& geomfields, 
		  PhononMacro& macro, Tkspace& kspace, PhononBCMap& bcMap,
		  const IntArray& BCArray, const BoolArray& BCfArray):
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
    _fgFinder()
    {}

  void COMETSolve()
  {
    const int cellcount=_cells.getSelfCount();
    
    for(int c=0;c<cellcount;c++)
      {	
	const int totalmodes=_kspace.gettotmodes();
	
	if(_BCArray[c]==0)  //not a reflecting boundary
	  {
	    TArrow AHMat(totalmodes+1);
	    TArray Bvec(totalmodes+1);
	    TArray Resid(totalmodes+1);
	    Bvec.zero();
	    
	    COMETConvection(c,AHMat,Bvec);
	    COMETCollision(c,AHMat,Bvec);
	    COMETEquilibrium(c,AHMat,Bvec);
	    Resid=Bvec;
	    AHMat.Solve(Bvec);
	    
	    Distribute(c,Bvec,Resid);
	  }
	else if(_BCArray[c]==1) //reflecting boundary
	  {
	    TSquare AMat(totalmodes+1);
	    TArray Bvec(totalmodes+1);
	    TArray Resid(totalmodes+1);
	    Bvec.zero();

	    throw CException("Haven't put it in yet.");
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");
      }   
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
	
	if(_BCfArray[f]==true) //If the face in question is a reflecting face
	  {
	    int Fgid=findFgId(f);
	    T refl=(*(_bcMap[Fgid]))["specifiedReflection"];

	    //first sweep - have to make sumVdotA
	    T sumVdotA=0.;
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    T dk3=mode.getdk3();
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
		      {//have to move through all other directions
			int ccount=1;
			for(int kk=0;kk<klen;kk++)
			  {
			    Tkvol& kkvol=_kspace.getkvol(kk);
			    const int numMODES=kkvol.getmodenum();
			    for(int mm=0;mm<numMODES;mm++)
			      {
				if(kk!=k && m!=mm)
				  {
				    VectorT3 vvg=kkvol.getmode(mm).getv();
				    T ddk3=kkvol.getmode(mm).getdk3();
				    T VdotA=vvg[0]*Af[0]+vvg[1]*Af[1]+vvg[2]*Af[2];    
				    if(VdotA>T_Scalar(0))
				      {
					Field& eefield=kvol.getmode(mm).getfield();
					TArray& eeArray=dynamic_cast<TArray&>(eefield[_cells]);
					Amat(count,ccount)-=flux*VdotA*ddk3;
					BVec[count-1]-=flux*VdotA*ddk3*eeArray[cell];
				      }
				  }
				ccount++;
			      }
			  }
		      }
		    count++;
		  }
	      }

	    //Second sweep - Now I have to divide some coeffs by sumVdotA
	    count=1;
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
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
			    for(int mm=0;mm<numMODES;mm++)
			      {
				if(kk!=k && m!=mm)
				  {
				    VectorT3 vvg=kkvol.getmode(mm).getv();
				    T VdotA=vvg[0]*Af[0]+vvg[1]*Af[1]+vvg[2]*Af[2];
				    if(VdotA>T_Scalar(0))
				      Amat(count,ccount)/=sumVdotA;
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
  
  void COMETCollision(const int cell, TArrow& Amat, TArray& BVec)
  {
    /* This is the COMET discretization for cells that do not
       do not have a face which is a reflecting boundary.  When
       there is a face with a reflecting boundary, we can no longer
       use an arrowhead matrix as the structure becomes unknown a priori.
     */
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    TArray& e0_Array=dynamic_cast<TArray&>(_macro.e0[_cells]);
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
	    Field& efield=mode.getfield();
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]);  
	    coeff=_cellVolume[cell]/tau;
	    Amat.getElement(count,order)+=coeff;
	    Amat.getElement(count,count)-=coeff;
	    BVec[count-1]-=coeff*eArray[cell];
	    BVec[count-1]+=coeff*e0_Array[cell];
	    count+=1;
	  }
      }
    
  }

  void COMETEquilibrium(const int cell, TArrow& Amat, TArray& BVec)
  {
    /* This is the COMET discretization for cells that do not
       do not have a face which is a reflecting boundary.  When
       there is a face with a reflecting boundary, we can no longer
       use an arrowhead matrix as the structure becomes unknown a priori.
     */
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    const T tauTot=_kspace.calcTauTot();
    TArray& e0_Array=dynamic_cast<TArray&>(_macro.e0[_cells]);
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
	    TArray& eArray=dynamic_cast<TArray&>(efield[_cells]); 
	    Amat.getElement(order,count)-=coeff;
	    BVec[totalmodes]-=coeff*eArray[cell];
	    count+=1;
	  }
      }
    Amat.getElement(order,order)=1.;
    BVec[totalmodes]+=e0_Array[cell];
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
    
    TArray& e0Array=dynamic_cast<TArray&>(_macro.e0[_cells]);
    e0Array[cell]-=BVec[totalmodes];
    TArray& e0ResArray=dynamic_cast<TArray&>(_macro.e0Residual[_cells]);
    e0ResArray[cell]=-Rvec[totalmodes];
  }

  void findResid()
  {
    const int cellcount=_cells.getSelfCount();
    const int totalmodes=_kspace.gettotmodes();
    TArray ResidSum(totalmodes+1);
    T ResidScalar=0.;
    ResidSum.zero();
    TArrow AHMat(totalmodes+1);

    for(int c=0;c<cellcount;c++)
      {	
	AHMat.zero();
	TArray Bvec(totalmodes+1);
	TArray Resid(totalmodes+1);
	Bvec.zero();
	Resid.zero();
	
	COMETConvection(c,AHMat,Bvec);
       	COMETCollision(c,AHMat,Bvec);
	COMETEquilibrium(c,AHMat,Bvec);
	Resid+=Bvec;
	Bvec.zero();
	ResidSum+=Resid;

	Distribute(c,Bvec,Resid);
      }

    for(int o=0;o<totalmodes;o++)
      {
	ResidSum[o]/=cellcount;
	ResidSum[o]/=totalmodes;
	ResidScalar+=ResidSum[o];
      }

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
  const TArray& _cellVolume;
  const VectorT3Array& _faceArea;
  PhononMacro& _macro;
  Tkspace& _kspace;
  PhononBCMap& _bcMap;
  const IntArray& _BCArray;
  const BoolArray& _BCfArray;
  T _aveResid;
  T _residChange;
  FaceToFg _fgFinder;
  
};


#endif
