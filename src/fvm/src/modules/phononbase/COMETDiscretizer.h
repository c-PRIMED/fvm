// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

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
#include "GradientModel.h"
#include "FluxLimiters.h"
#include <omp.h>

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
  typedef KSConnectivity<T> TKConnectivity;
  typedef TKConnectivity* TKCptr;
  typedef vector<TKCptr> TKClist;
  typedef map<int, TKClist> FgTKClistMap;
  typedef Gradient<T> GradType;
  typedef Array<GradType> GradArray;
  typedef GradientModel<T> GradModelType;
  typedef typename GradModelType::GradMatrixType GradMatrix;

 COMETDiscretizer(const Mesh& mesh, const GeomFields& geomfields, 
		  PhononMacro& macro, Tkspace& kspace, COMETBCMap& bcMap,
		  const IntArray& BCArray, const IntArray& BCfArray, COpts& options,
		  const FgTKClistMap& FgToKsc):
  _mesh(mesh),
    _geomFields(geomfields),
    _cells(mesh.getCells()),
    _faces(mesh.getFaces()),
    _cellFaces(mesh.getCellFaces()),
    _faceCells(mesh.getAllFaceCells()),
    _faceAreaMag((dynamic_cast<const TArray&>(_geomFields.areaMag[_faces]))),
    _faceArea(dynamic_cast<const VectorT3Array&>(_geomFields.area[_faces])),
    _cellVolume(dynamic_cast<const TArray&>(_geomFields.volume[_cells])),
    _cellCoords(dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_cells])),
    _macro(macro),
    _kspace(kspace),
    _bcMap(bcMap),
    _BCArray(BCArray),
    _BCfArray(BCfArray),
    _aveResid(-1.),
    _residChange(-1.),
    _fgFinder(),
    _options(options),
    _FaceToKSC(FgToKsc),
    _eArray(kspace.geteArray()),
    _e0Array(kspace.gete0Array()),
    _resArray(kspace.getResArray())
    {}

  void COMETSolveFine(const int dir,const int level)
  {
    const int cellcount=_cells.getSelfCount();
    int start;

    if(dir==1)
      start=0;
    if(dir==-1)
      start=cellcount-1;
    const int totalmodes=_kspace.gettotmodes();
    TArray Bvec(totalmodes+1);
    TArray Resid(totalmodes+1);
    TArrow AMat(totalmodes+1);

    const GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);
    
    for(int c=start;((c<cellcount)&&(c>-1));c+=dir)
      {	
		
	if(_BCArray[c]==0)  //no reflections at all--interior cell or temperature boundary
	  {
	    T dt=1;
	    int NewtonIters=0;
	    //updateGhostFine(c, gradMatrix);
	  
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {
		
		Bvec.zero();
		Resid.zero();
		AMat.zero();
		
		COMETConvectionFine(c,AMat,Bvec,gradMatrix);
		COMETCollision(c,&AMat,Bvec);
		COMETEquilibrium(c,&AMat,Bvec);

		if(_options.withNormal)
		  COMETShifted(c,&AMat,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
		
		Resid=Bvec;
		AMat.Solve(Bvec);

		Distribute(c,Bvec,Resid);
		updatee0(c);
		updateGhostFine(c, gradMatrix);
		NewtonIters++;
		dt=fabs(Bvec[totalmodes]);
	      }
	  }
	else if(_BCArray[c]==1) //Implicit reflecting boundary only
	  {
	    T dt=1;
	    int NewtonIters=0;
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {
		TSquare AMatS(totalmodes+1);
		Bvec.zero();
		Resid.zero();
		AMatS.zero();
		
		COMETConvection(c,AMatS,Bvec);
		COMETCollision(c,&AMatS,Bvec);
		COMETEquilibrium(c,&AMatS,Bvec);
		  
		if(level>0)
		  addFAS(c,Bvec);
	       
		Resid=Bvec;
		AMatS.Solve(Bvec);
		
		Distribute(c,Bvec,Resid);
		updatee0(c);
		NewtonIters++;
		dt=fabs(Bvec[totalmodes]);
	      }
	  }
	else if(_BCArray[c]==2)  //Explicit boundary only
	  {
	    T dt=1;
	    int NewtonIters=0;
	    updateGhostFine(c, gradMatrix);
	    //updateGhostCoarse(c);
	  
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {

		Bvec.zero();
		Resid.zero();
		AMat.zero();
		
		COMETConvectionFine(c,AMat,Bvec,gradMatrix);
		COMETCollision(c,&AMat,Bvec);
		COMETEquilibrium(c,&AMat,Bvec);
		
		if(_options.withNormal)
		  COMETShifted(c,&AMat,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
		
		Resid=Bvec;
		AMat.Solve(Bvec);

		Distribute(c,Bvec,Resid);
		dt=fabs(Bvec[totalmodes]);
		updatee0(c);
		updateGhostFine(c, gradMatrix);
		if(!_FaceToKSC.empty())
		  correctInterface(c,Bvec);
		NewtonIters++;
	      }
	  }
	else if(_BCArray[c]==3) //Mix Implicit/Explicit reflecting boundary
	  {
	    T dt=1;
	    int NewtonIters=0;
	    updateGhostFine(c,gradMatrix);
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {
		TSquare AMatS(totalmodes+1);
		Bvec.zero();
		Resid.zero();
		AMatS.zero();
		
		COMETConvection(c,AMatS,Bvec);
		COMETCollision(c,&AMatS,Bvec);
		COMETEquilibrium(c,&AMatS,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
	       
		Resid=Bvec;
		AMatS.Solve(Bvec);
		
		Distribute(c,Bvec,Resid);
		updatee0(c);
		updateGhostFine(c, gradMatrix);
		NewtonIters++;
		dt=fabs(Bvec[totalmodes]);
	      }
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");

      }
    
  }

  void COMETSolveCoarse(const int dir,const int level)
  {
    const int cellcount=_cells.getSelfCount();
    int start;

    if(dir==1)
      start=0;
    if(dir==-1)
      start=cellcount-1;
    const int totalmodes=_kspace.gettotmodes();
    TArray Bvec(totalmodes+1);
    TArray Resid(totalmodes+1);
    TArrow AMat(totalmodes+1);
    
    for(int c=start;((c<cellcount)&&(c>-1));c+=dir)
      {	
		
	if(_BCArray[c]==0)  //no reflections at all--interior cell or temperature boundary
	  {
	    T dt=1;
	    int NewtonIters=0;
	  
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {
		
		Bvec.zero();
		Resid.zero();
		AMat.zero();
		
		COMETConvectionCoarse(c,AMat,Bvec);
		COMETCollision(c,&AMat,Bvec);
		COMETEquilibrium(c,&AMat,Bvec);

		if(_options.withNormal)
		  COMETShifted(c,&AMat,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
		
		Resid=Bvec;
		AMat.Solve(Bvec);

		Distribute(c,Bvec,Resid);
		updatee0(c);
		NewtonIters++;
		dt=fabs(Bvec[totalmodes]);
	      }
	  }
	else if(_BCArray[c]==1) //Implicit reflecting boundary only
	  {
	    T dt=1;
	    int NewtonIters=0;
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {
		TSquare AMatS(totalmodes+1);
		Bvec.zero();
		Resid.zero();
		AMatS.zero();
		
		COMETConvection(c,AMatS,Bvec);
		COMETCollision(c,&AMatS,Bvec);
		COMETEquilibrium(c,&AMatS,Bvec);
		  
		if(level>0)
		  addFAS(c,Bvec);
	       
		Resid=Bvec;
		AMatS.Solve(Bvec);
		
		Distribute(c,Bvec,Resid);
		updatee0(c);
		NewtonIters++;
		dt=fabs(Bvec[totalmodes]);
	      }
	  }
	else if(_BCArray[c]==2)  //Explicit boundary only
	  {
	    T dt=1;
	    int NewtonIters=0;
	    updateGhostCoarse(c);
	  
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {

		Bvec.zero();
		Resid.zero();
		AMat.zero();
		
		COMETConvectionCoarse(c,AMat,Bvec);
		COMETCollision(c,&AMat,Bvec);
		COMETEquilibrium(c,&AMat,Bvec);
		
		if(_options.withNormal)
		  COMETShifted(c,&AMat,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
		
		Resid=Bvec;
		AMat.Solve(Bvec);

		Distribute(c,Bvec,Resid);
		dt=fabs(Bvec[totalmodes]);
		updatee0(c);
		updateGhostCoarse(c);
		if(!_FaceToKSC.empty())
		  correctInterface(c,Bvec);
		NewtonIters++;
	      }
	  }
	else if(_BCArray[c]==3) //Mix Implicit/Explicit reflecting boundary
	  {
	    T dt=1;
	    int NewtonIters=0;
	    updateGhostCoarse(c);
	    while(dt>_options.NewtonTol && NewtonIters<50)
	      {
		TSquare AMatS(totalmodes+1);
		Bvec.zero();
		Resid.zero();
		AMatS.zero();
		
		COMETConvection(c,AMatS,Bvec);
		COMETCollision(c,&AMatS,Bvec);
		COMETEquilibrium(c,&AMatS,Bvec);
		
		if(level>0)
		  addFAS(c,Bvec);
	       
		Resid=Bvec;
		AMatS.Solve(Bvec);
		
		Distribute(c,Bvec,Resid);
		updatee0(c);
		updateGhostCoarse(c);
		NewtonIters++;
		dt=fabs(Bvec[totalmodes]);
	      }
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");

      }
    
  }

  void COMETConvectionFine(const int cell0, TArrow& Amat, 
			   TArray& BVec, const GradMatrix& gMat)
  {

    const int neibcount=_cellFaces.getCount(cell0);
    GradArray Grads(_kspace.gettotmodes());
    TArray pointMin(_kspace.gettotmodes());
    pointMin=-1;
    TArray pointMax(_kspace.gettotmodes());
    pointMax.zero();
    TArray neibMin(_kspace.gettotmodes());
    neibMin=-1;
    TArray neibMax(_kspace.gettotmodes());
    neibMax.zero();
    TArray pointLim(_kspace.gettotmodes());
    pointLim=100.;
    TArray neibLim(_kspace.gettotmodes());
    neibLim=100.;

    for(int k=0;k<_kspace.gettotmodes();k++)
      Grads[k].zero();

    const VectorT3Array& faceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_faces]);

    for(int j=0;j<neibcount;j++)    //first loop to get grad and max/min vals
      {
	const int f=_cellFaces(cell0,j);
	int cell1=_faceCells(f,1);
	if(cell1==cell0)
	  cell1=_faceCells(f,0);
	const VectorT3 Gcoeff=gMat.getCoeff(cell0, cell1);

	int c0ind=_kspace.getGlobalIndex(cell0,0);
	int c1ind=_kspace.getGlobalIndex(cell1,0);
	
	for(int k=0;k<_kspace.gettotmodes();k++)
	  {
	    const T e1=_eArray[c1ind];
	    const T e0=_eArray[c0ind];
	    T& curMax=pointMax[k];
	    T& curMin=pointMin[k];
	    Grads[k].accumulate(Gcoeff, e1-e0);

	    if(e1>curMax)
	      curMax=e1;
	    
	    if(curMin==-1)
	      curMin=e1;
	    else
	      {
		if(e1<curMin)
		  curMin=e1;
	      }

	    c0ind++;
	    c1ind++;
	  }  
      }

    for(int j=0;j<neibcount;j++)    //second loop to calculate limiting coeff
      {
	const int f=_cellFaces(cell0,j);
	int cell1=_faceCells(f,1);
	if(cell1==cell0)
	  cell1=_faceCells(f,0);

	const VectorT3 dr0(faceCoords[f]-_cellCoords[cell0]);
	int c0ind=_kspace.getGlobalIndex(cell0,0);
	//vanLeer lf;
	
	for(int k=0;k<_kspace.gettotmodes();k++)
	  {
	    const T minVal=pointMin[k];
	    const T maxVal=pointMax[k];
	    const T de0(Grads[k]*dr0);
	    T& cl=pointLim[k];
	    //computeLimitCoeff(cl, _eArray[c0ind], de0, minVal, maxVal, lf);
	    c0ind++;
	  }  
      }

    for(int j=0;j<neibcount;j++)  //loop to construct point matrix
      {
	const int f=_cellFaces(cell0,j);
	int cell1=_faceCells(f,1);
	VectorT3 Af=_faceArea[f];
	if(cell1==cell0)
	  {
	    cell1=_faceCells(f,0);
	    Af*=-1.;
	  }
	const int klen=_kspace.getlength();

	GradArray NeibGrads(_kspace.gettotmodes());

	for(int k=0;k<_kspace.gettotmodes();k++)
	  NeibGrads[k].zero();
	
	neibMax.zero();
	neibMin=-1.;
	const int neibcount1=_cellFaces.getCount(cell1);
	for(int nj=0;nj<neibcount1;nj++)  //loop to make gradients and min/max
	  {
	    const int f1=_cellFaces(cell1,nj);
	    int cell11=_faceCells(f1,1);
	    if(cell1==cell11)
	      cell11=_faceCells(f1,0);
	
	    const VectorT3 Gcoeff=gMat.getCoeff(cell1, cell11);

	    int c1ind=_kspace.getGlobalIndex(cell1,0);
	    int c11ind=_kspace.getGlobalIndex(cell11,0);
	
	    for(int k=0;k<_kspace.gettotmodes();k++)
	      {
		const T e1=_eArray[c1ind];
		const T e11=_eArray[c11ind];
		T& curMax=neibMax[k];
		T& curMin=neibMin[k];
		NeibGrads[k].accumulate(Gcoeff,e11-e1);

		if(e11>curMax)
		  curMax=e11;
	    
		if(curMin==-1)
		  curMin=e11;
		else
		  {
		    if(e11<curMin)
		      curMin=e11;
		  }

		c11ind++;
		c1ind++;
	      }
	  }
	  
	neibLim=100.;
	for(int nj=0;nj<neibcount1;nj++)  //second loop to calculate limiting coeff
	  {
	    const int f1=_cellFaces(cell1,nj);
	    int cell11=_faceCells(f1,1);
	    if(cell1==cell11)
	      cell11=_faceCells(f1,0);
	
	    const VectorT3 dr1(faceCoords[f1]-_cellCoords[cell1]);
	    int c1ind=_kspace.getGlobalIndex(cell1,0);
	    //vanLeer lf;
	
	    for(int k=0;k<_kspace.gettotmodes();k++)
	      {
		const T minVal=neibMin[k];
		const T maxVal=neibMax[k];
		const T de1(NeibGrads[k]*dr1);
		T& cl=neibLim[k];
		//computeLimitCoeff(cl, _eArray[c1ind], de1, minVal, maxVal, lf);
		if(_BCfArray[f]!=0)
		  NeibGrads[k].zero();
		c1ind++;
	      }
	  }


	T flux;

	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		const int count=mode.getIndex();
		VectorT3 vg=mode.getv();
		T tau=mode.gettau();
		flux=vg[0]*Af[0]+vg[1]*Af[1]+vg[2]*Af[2];

		const int c0ind=_kspace.getGlobalIndex(cell0,count-1);
		const int c1ind=_kspace.getGlobalIndex(cell1,count-1);

		const GradType& grad=Grads[count-1];
		const GradType& neibGrad=NeibGrads[count-1];
		
		if(flux>T_Scalar(0))
		  {
		    const VectorT3 rVec=_cellCoords[cell1]-_cellCoords[cell0];
		    const VectorT3 fVec=faceCoords[f]-_cellCoords[cell0];
		    
		    T SOU=(fVec[0]*grad[0]+fVec[1]*grad[1]+
			   fVec[2]*grad[2])*pointLim[count-1];
		    Amat.getElement(count,count)-=flux;
		    BVec[count-1]-=flux*_eArray[c0ind]-flux*SOU;
		  }
		else
		  {
		    const VectorT3 rVec=_cellCoords[cell0]-_cellCoords[cell1];
		    const VectorT3 fVec=faceCoords[f]-_cellCoords[cell1];
		    
		    T SOU=(fVec[0]*neibGrad[0]+fVec[1]*neibGrad[1]+
			   fVec[2]*neibGrad[2])*neibLim[count-1];
		    BVec[count-1]-=flux*_eArray[c1ind]-flux*SOU;
		  }
		
	      }
	  }
      }

    /*
    const int klen=_kspace.getlength();
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex();
	    VectorT3 vg=mode.getv();
	    T tau=mode.gettau();
	    T vmag=pow(vg.mag2(),0.5);
	    T mfp=vmag*tau;
	    BVec[count-1]*=mfp;
	    Amat.getElement(count,count)*=mfp;
	  }
      }*/
    
  }

  void COMETConvectionCoarse(const int cell0, TArrow& Amat, TArray& BVec)
  {
    /* This is the COMET discretization for cells that do not
       do not have a face which is a reflecting boundary.  When
       there is a face with a reflecting boundary, we can no longer
       use an arrowhead matrix as the structure becomes unknown a priori.
    */
    const int neibcount=_cellFaces.getCount(cell0);
    
    for(int j=0;j<neibcount;j++)
      {
	const int f=_cellFaces(cell0,j);
	int cell1=_faceCells(f,1);
	VectorT3 Af=_faceArea[f];
	if(cell1==cell0)
	  {
	    cell1=_faceCells(f,0);
	    Af*=-1.;
	  }
	const int klen=_kspace.getlength();
	
	T flux;
	T fluxX(0);
	T fluxY(0);

	for(int k=0;k<klen;k++)
	  {

	    Tkvol& kvol=_kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		const int count=mode.getIndex();
		const VectorT3 vg=mode.getv();
		const T vmag=pow(vg.mag2(),0.5);
		const T tau=mode.gettau();
		const T mfp=1.;//vmag*tau;
		flux=vg[0]*(Af[0]/mfp)+vg[1]*(Af[1]/mfp)+vg[2]*(Af[2]/mfp);
		fluxX+=vg[0]*(Af[0]/mfp);
		fluxY+=vg[1]*(Af[1]/mfp);
		
		if(flux>T_Scalar(0))
		  {
		    Amat.getElement(count,count)-=flux;
		    BVec[count-1]-=flux*_eArray[_kspace.getGlobalIndex(cell0,count-1)];
		  }
		else
		  BVec[count-1]-=flux*_eArray[_kspace.getGlobalIndex(cell1,count-1)];
	      }
	  }

	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		const int count=mode.getIndex();
		VectorT3 vg=mode.getv();
		T vmag=pow(vg.mag2(),0.5);
		T tau=mode.gettau();
		T mfp=1.;//vmag*tau;
		BVec[count-1]*=mfp;
		Amat.getElement(count,count)*=mfp;
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

	if(cell2==cell)
	  {
	    Af=Af*(-1.);
	    cell2=_faceCells(f,0);
	  }
       
	if(_BCfArray[f]==2) //If the face in question is an implicit reflecting face
	  {
	    int Fgid=findFgId(f);
	    T refl=(*(_bcMap[Fgid]))["specifiedReflection"];
	    T oneMinusRefl=1.-refl;

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
		    const int count=mode.getIndex();
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
		      {
			sumVdotA-=flux*dk3;
			Refl_pair& rpairs=mode.getReflpair(Fgid);
			const int kk=rpairs.second.second;
			Tkvol& kkvol=_kspace.getkvol(kk);
			const int ccount=kkvol.getmode(m).getIndex();
			Field& eefield=kkvol.getmode(m).getfield();
			TArray& eeArray=dynamic_cast<TArray&>(eefield[_cells]);
			Amat(count,ccount)-=flux*refl;
			BVec[count-1]-=eeArray[cell]*refl*flux;
		      }
		  }
	      }
	    
	    T inverseSumVdotA=1./sumVdotA;

	    //Second sweep
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    const int count=kvol.getmode(m).getIndex();
		    VectorT3 vg=kvol.getmode(m).getv();
		    flux=vg[0]*Af[0]+vg[1]*Af[1]+vg[2]*Af[2];
		    if(flux<T_Scalar(0))
		      {
			for(int kk=0;kk<klen;kk++)
			  {
			    Tkvol& kkvol=_kspace.getkvol(kk);
			    const int numMODES=kkvol.getmodenum();
			    T ddk3=kkvol.getdk3();
			    for(int mm=0;mm<numMODES;mm++)
			      {
				VectorT3 vvg=kkvol.getmode(mm).getv();
				T VdotA=vvg[0]*Af[0]+vvg[1]*Af[1]+vvg[2]*Af[2];
				if(VdotA>T_Scalar(0))
				  {
				    const int ccount=kkvol.getmode(mm).getIndex();
				    Field& eefield=kkvol.getmode(mm).getfield();
				    TArray& eeArray=dynamic_cast<TArray&>(eefield[_cells]);
				    Amat(count,ccount)-=flux*VdotA*ddk3*oneMinusRefl*inverseSumVdotA;
				    BVec[count-1]-=flux*VdotA*ddk3
				      *eeArray[cell]*inverseSumVdotA*oneMinusRefl;
				  }
			      }
			  }
		      }
		  }
	      }
	    
	    
	  }
	else  //if the face in question is not implicit
	  {
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    const int count=mode.getIndex();
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
    T coeff;
    int cellIndex=_kspace.getGlobalIndex(cell,0);
    
#pragma omp parallel for shared(BVec, Amat, _kspace) private(k, m, numModes, count, mode, kvol, tau, de0dT)
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex();
	    T tau=mode.gettau();
	    T de0dT=mode.calcde0dT(Tlold[cell]);
	    coeff=_cellVolume[cell]/tau;
	    Amat->getElement(count,order)+=coeff*de0dT;
	    Amat->getElement(count,count)-=coeff;
	    BVec[count-1]-=coeff*_eArray[cellIndex];
	    BVec[count-1]+=coeff*_e0Array[cellIndex];
	    cellIndex++;
	  }
      }
    
  }

  void COMETEquilibrium(const int cell, TMatrix* Amat, TArray& BVec)
  {
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    TArray& Tlold=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    const T tauTot=_kspace.getde0taudT(Tlold[cell]);//getde0taudTgray();
    T coeff;
    const T DK3=_kspace.getDK3();
    int cellIndex=_kspace.getGlobalIndex(cell,0);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const T dk3=kvol.getdk3();
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex();
	    T tau=mode.gettau();
	    coeff=(dk3/DK3)/tau/tauTot;
	    Amat->getElement(order,count)+=coeff;
	    BVec[totalmodes]+=coeff*_eArray[cellIndex];
	    BVec[totalmodes]-=coeff*_e0Array[cellIndex];
	    cellIndex++;
	  }
      }
    Amat->getElement(order,order)=-1.;
    BVec[totalmodes]*=DK3;
  }

  void COMETShifted(const int cell, TMatrix* Amat, TArray& BVec)
  { //adds to collision and equilibrium
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    const int order=totalmodes+1;
    TArray& Tlold=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    const T tauTot=_kspace.getde0taudT(Tlold[cell]);
    T coeff;
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex();
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
	  }
      }
    
  }

  void Distribute(const int cell, TArray& BVec, TArray& Rvec)
  {
    const int klen=_kspace.getlength();
    const int totalmodes=_kspace.gettotmodes();
    int cellIndex=_kspace.getGlobalIndex(cell,0);

    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex()-1;
	    _eArray[cellIndex]-=BVec[count];
	    _resArray[cellIndex]=-Rvec[count];
	    cellIndex++;
	  }
      }
    
    TArray& TlArray=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    TlArray[cell]-=BVec[totalmodes];
    TArray& deltaTArray=dynamic_cast<TArray&>(_macro.deltaT[_cells]);
    deltaTArray[cell]=BVec[totalmodes];
    TArray& TlResArray=dynamic_cast<TArray&>(_macro.TlResidual[_cells]);
    TlResArray[cell]=-Rvec[totalmodes];
  }

  void findResidFine()
  {
    const int cellcount=_cells.getSelfCount();
    const int totalmodes=_kspace.gettotmodes();
    TArray ResidSum(totalmodes+1);
    TArray Bsum(totalmodes+1);
    T ResidScalar=0.;
    T traceSum=0.;
    TArray Bvec(totalmodes+1);
    TArray Resid(totalmodes+1);
    TArray Dummy(totalmodes+1);
    TArrow AMat(totalmodes+1);
    ResidSum.zero();
    Bsum.zero();
    Dummy.zero();

    GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);
    
    for(int c=0;c<cellcount;c++)
      {	
	Bvec.zero();
	Resid.zero();
	Dummy.zero();

	if(_BCArray[c]==0 || _BCArray[c]==2)  //Arrowhead
	  {
	      
	    updateGhostFine(c, gradMatrix);
	    //updateGhostCoarse(c);

	    AMat.zero();

	    COMETConvectionFine(c,AMat,Bvec,gradMatrix);	
	    COMETCollision(c,&AMat,Bvec);
	    COMETEquilibrium(c,&AMat,Bvec);
	    
	    traceSum+=AMat.getTraceAbs();
	    
	    Resid=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);
	    ArrayAbs(Resid);
	    ResidSum+=Resid;

	    AMat.multiply(Resid,Bvec);
	    Resid=Bvec;
	    
	    makeValueArray(c,Bvec);
	    AMat.multiply(Bvec,Dummy);
	    Bvec=Dummy;
	    ArrayAbs(Bvec);
	    ArrayAbs(Resid);
	    Bsum+=Bvec;	    
	  }
	else if(_BCArray[c]==1 || _BCArray[c]==3) //General Dense
	  {
	    TSquare AMatS(totalmodes+1);
	    AMatS.zero();
	    
	    COMETConvection(c,AMatS,Bvec);	
	    COMETCollision(c,&AMatS,Bvec);
	    COMETEquilibrium(c,&AMatS,Bvec);
	    
	    traceSum+=AMatS.getTraceAbs();
	    Resid=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);

	    AMatS.multiply(Resid,Bvec);
	    Resid=Bvec;

	    makeValueArray(c,Bvec);
	    ArrayAbs(Resid);
	    ArrayAbs(Bvec);
	    Bsum+=Bvec;
	    ResidSum+=Resid;
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");
      }

    //traceSum=0;  //added
    for(int o=0;o<totalmodes+1;o++)
      {
	ResidScalar+=ResidSum[o];
	//traceSum+=Bsum[o]; //added
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

  void findResidCoarse(const bool plusFAS)
  {
    const int cellcount=_cells.getSelfCount();
    const int totalmodes=_kspace.gettotmodes();
    TArray ResidSum(totalmodes+1);
    TArray Bsum(totalmodes+1);
    T ResidScalar=0.;
    T traceSum=0.;
    TArray Bvec(totalmodes+1);
    TArray Resid(totalmodes+1);
    TArray Dummy(totalmodes+1);
    TArrow AMat(totalmodes+1);
    ResidSum.zero();
    Bsum.zero();
    Dummy.zero();
    
    for(int c=0;c<cellcount;c++)
      {	
	Bvec.zero();
	Resid.zero();
	Dummy.zero();

	if(_BCArray[c]==0 || _BCArray[c]==2)  //Arrowhead
	  {
	    if(_BCArray[c]==2)
	      updateGhostCoarse(c);

	    AMat.zero();

	    COMETConvectionCoarse(c,AMat,Bvec);	
	    COMETCollision(c,&AMat,Bvec);
	    COMETEquilibrium(c,&AMat,Bvec);
	    
	    if(plusFAS)
	      addFAS(c,Bvec);

	    traceSum+=AMat.getTraceAbs();
	    
	    Resid=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);
	    ArrayAbs(Resid);
	    ResidSum+=Resid;

	    AMat.multiply(Resid,Bvec);
	    Resid=Bvec;
	    
	    makeValueArray(c,Bvec);
	    AMat.multiply(Bvec,Dummy);
	    Bvec=Dummy;
	    ArrayAbs(Bvec);
	    ArrayAbs(Resid);
	    Bsum+=Bvec;	    
	  }
	else if(_BCArray[c]==1 || _BCArray[c]==3) //General Dense
	  {
	    TSquare AMatS(totalmodes+1);
	    AMatS.zero();
	    
	    COMETConvection(c,AMatS,Bvec);	
	    COMETCollision(c,&AMatS,Bvec);
	    COMETEquilibrium(c,&AMatS,Bvec);
	    
	    if(plusFAS)
	      addFAS(c,Bvec);

	    traceSum+=AMatS.getTraceAbs();
	    Resid=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);

	    AMatS.multiply(Resid,Bvec);
	    Resid=Bvec;

	    makeValueArray(c,Bvec);
	    ArrayAbs(Resid);
	    ArrayAbs(Bvec);
	    Bsum+=Bvec;
	    ResidSum+=Resid;
	  }
	else
	  throw CException("Unexpected value for boundary cell map.");
      }

    //traceSum=0;  //added
    for(int o=0;o<totalmodes+1;o++)
      {
	ResidScalar+=ResidSum[o];
	//traceSum+=Bsum[o]; //added
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

    TArray rArray(order);
    
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex();
	    Field& resfield=mode.getresid();
	    TArray& resArray=dynamic_cast<TArray&>(resfield[_cells]);
	    rArray[count]=resArray[c];
	  }
      }
    rArray[order]=_macro.TlResidual[c];
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
    cout<<"Face index: "<<faceIndex<<endl;
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
    const int totmodes=_kspace.gettotmodes();
    int cellIndex=_kspace.getGlobalIndex(c,0);
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    const int count=mode.getIndex()-1;
	    o[count]=_eArray[cellIndex];
	    cellIndex++;
	  }
      }
    TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    o[totmodes]=Tl[c];
  }

  void addFAS(const int c, TArray& bVec)
  {
    const int totmodes=_kspace.gettotmodes();
    _kspace.addFAS(c,bVec);
    TArray& fasArray=dynamic_cast<TArray&>(_macro.TlFASCorrection[_cells]);
    bVec[totmodes]-=fasArray[c];
  }

  void updatee0()
  {
    TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);
   
    const T cellCount=_cells.getSelfCount();
    int klen=_kspace.getlength();
    for(int c=0;c<cellCount;c++)
      {
	int cellIndex=_kspace.getGlobalIndex(c,0);
	for(int k=0;k<klen;k++)
	  {
	    Tkvol& kvol=_kspace.getkvol(k);
	    const int numModes=kvol.getmodenum();
	    for(int m=0;m<numModes;m++)
	      {
		Tmode& mode=kvol.getmode(m);
		_e0Array[cellIndex]=mode.calce0(Tl[c]);
		cellIndex++;
	      }
	  }
      }
  }

  void updatee0(const int c)
  {
    TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);

    int klen=_kspace.getlength();
    int cellIndex=_kspace.getGlobalIndex(c,0);
    for(int k=0;k<klen;k++)
      {
	Tkvol& kvol=_kspace.getkvol(k);
	const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    _e0Array[cellIndex]=mode.calce0(Tl[c]);
	    cellIndex++;
	  }
      }
  }

  void updateGhostFine(const int cell, const GradMatrix& gMat)
  {
    const int neibcount=_cellFaces.getCount(cell);
    const int klen=_kspace.getlength();
    //TArray SumEVdotA(neibcount);
    //TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    //T DK3=_kspace.getDK3();

    //SumEVdotA.zero();

    const VectorT3Array& faceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_faces]);

    TArray pointMin(_kspace.gettotmodes());
    pointMin=-1;
    TArray pointMax(_kspace.gettotmodes());
    pointMax.zero();
    TArray pointLim(_kspace.gettotmodes());
    pointLim=100.;

    GradArray Grads(_kspace.gettotmodes());
    Grads.zero();

    VectorT3 Gcoeff;

    for(int j=0;j<neibcount;j++)    //first loop to get grad and max/min vals
      {
	const int f=_cellFaces(cell,j);
	int cell1=_faceCells(f,1);
	if(cell1==cell)
	  cell1=_faceCells(f,0);
	const VectorT3 Gcoeff=gMat.getCoeff(cell, cell1);

	int c0ind=_kspace.getGlobalIndex(cell,0);
	int c1ind=_kspace.getGlobalIndex(cell1,0);
	
	for(int k=0;k<_kspace.gettotmodes();k++)
	  {
	    const T e1=_eArray[c1ind];
	    const T e0=_eArray[c0ind];
	    T& curMax=pointMax[k];
	    T& curMin=pointMin[k];
	    Grads[k].accumulate(Gcoeff, e1-e0);

	    if(e1>curMax)
	      curMax=e1;
	    
	    if(curMin==-1)
	      curMin=e1;
	    else
	      {
		if(e1<curMin)
		  curMin=e1;
	      }

	    c0ind++;
	    c1ind++;
	  }  
      }

    for(int j=0;j<neibcount;j++)    //second loop to calculate limiting coeff
      {
	const int f=_cellFaces(cell,j);
	int cell1=_faceCells(f,1);
	if(cell1==cell)
	  cell1=_faceCells(f,0);

	const VectorT3 dr0(faceCoords[f]-_cellCoords[cell]);
	int c0ind=_kspace.getGlobalIndex(cell,0);
	//vanLeer lf;
	
	for(int k=0;k<_kspace.gettotmodes();k++)
	  {
	    const T minVal=pointMin[k];
	    const T maxVal=pointMax[k];
	    const T de0(Grads[k]*dr0);
	    T& cl=pointLim[k];
	    //computeLimitCoeff(cl, _eArray[c0ind], de0, minVal, maxVal, lf);
	    c0ind++;
	  }  
      }
    

    for(int j=0;j<neibcount;j++)
      {

	const int f=_cellFaces(cell,j);
	
	if(_BCfArray[f]==1)
	  {//temperature
	    int Fgid=findFgId(f);
	    int cell2=_faceCells(f,1);
	    VectorT3 Af=_faceArea[f];
		    
	    if(cell2==cell)
	      {
		Af=Af*-1.;
		cell2=_faceCells(f,0);
	      }

	    int cellIndex=_kspace.getGlobalIndex(cell,0);
	    int cell2Index=_kspace.getGlobalIndex(cell2,0);

	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		//T dk3=kvol.getdk3();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    const int count=mode.getIndex();
		    
		    GradType& grad=Grads[count-1];

		    if(VdotA>0)
		      {
			VectorT3 rVec=_cellCoords[cell2]-_cellCoords[cell];
			VectorT3 fVec=faceCoords[f]-_cellCoords[cell];
			T SOU=(fVec[0]*grad[0]+fVec[1]*grad[1]+
			   fVec[2]*grad[2])*pointLim[count-1];
			_eArray[cell2Index]=_eArray[cellIndex]+SOU;
		      }
		    cellIndex++;
		    cell2Index++;
		  }
	      }
      
	  }
	else if(_BCfArray[f]==3)
	  {//reflecting
	    int Fgid=findFgId(f);
	    const T refl=(*(_bcMap[Fgid]))["specifiedReflection"];
	    int cell2=_faceCells(f,1);
	    VectorT3 Af=_faceArea[f];
		    
	    if(cell2==cell)
	      {
		Af=Af*-1.;
		cell2=_faceCells(f,0);
	      }

	    const int cellstart=_kspace.getGlobalIndex(cell,0);
	    const int cell2start=_kspace.getGlobalIndex(cell2,0);
	    int cellIndex=_kspace.getGlobalIndex(cell,0);
	    int cell2Index=_kspace.getGlobalIndex(cell2,0);

	    //first loop to sum up the diffuse reflection
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		//T dk3=kvol.getdk3();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    const int count=mode.getIndex();
		    
		    GradType& grad=Grads[count-1];

		    if(VdotA>0)
		      {
			VectorT3 rVec=_cellCoords[cell2]-_cellCoords[cell];
			VectorT3 fVec=faceCoords[f]-_cellCoords[cell];
			T r=gMat.computeR(grad,_eArray,rVec,cellIndex,cell2Index);
			T SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);//*vanLeer(r);
			_eArray[cell2Index]=(_eArray[cellIndex]+SOU)*refl;
		      }
		    cellIndex++;
		    cell2Index++;
		  }
	      }

	    cellIndex=_kspace.getGlobalIndex(cell,0);
	    cell2Index=_kspace.getGlobalIndex(cell2,0);

	    //first loop to sum up the diffuse reflection
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		//T dk3=kvol.getdk3();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    const int count=mode.getIndex();
		    
		    if(VdotA<0)
		      {
			Refl_pair& rpairs=mode.getReflpair(Fgid);
			const int kk=rpairs.second.second;
			Tkvol& kkvol=_kspace.getkvol(kk);
			Tmode& refMode=kkvol.getmode(m);
			const int rIndex=cell2start+refMode.getIndex()-1;
			_eArray[cell2Index]=_eArray[rIndex]*refl;
			}
		    cellIndex++;
		    cell2Index++;
		  }
	      }

	  }
	else if(_BCfArray[f]==1)
	  {
	    int cell2=_faceCells(f,1);
	    VectorT3 Af=_faceArea[f];

	    if(cell2==cell)
	      {
		Af=Af*-1.;
		cell2=_faceCells(f,0);
	      }

	    const int cellstart=_kspace.getGlobalIndex(cell,0);
	    const int cell2start=_kspace.getGlobalIndex(cell2,0);
	    int cellIndex=_kspace.getGlobalIndex(cell,0);
	    int cell2Index=_kspace.getGlobalIndex(cell2,0);

	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    const int count=mode.getIndex();
		    
		    GradType& grad=Grads[count-1];

		    if(VdotA>0)
		      {
			VectorT3 rVec=_cellCoords[cell2]-_cellCoords[cell];
			VectorT3 fVec=faceCoords[f]-_cellCoords[cell];
			T r=gMat.computeR(grad,_eArray,rVec,cellIndex,cell2Index);
			T SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);//*vanLeer(r);
			_eArray[cell2Index]=(_eArray[cellIndex]+SOU);
		      }
		    cellIndex++;
		    cell2Index++;
		  }
	      }

	  }
      }

    /*
    T wallTemp(Tl[cell]);
    for(int j=0;j<neibcount;j++)
      {
	const int f=_cellFaces(cell,j);
	if(_BCfArray[f]==3)
	  {
	    int cell2=_faceCells(f,1);
	    VectorT3 Af=_faceArea[f];
	    if(cell2==cell)
	      {
		Af=Af*-1.;
		cell2=_faceCells(f,0);
	      }
	    T esum=SumEVdotA[j]*DK3;
	    _kspace.calcTemp(wallTemp,esum,Af);
	    SumEVdotA[j]=wallTemp;
	    Tl[cell2]=wallTemp;
	  }
      }

    //second loop to add in the diffuse reflection
    for(int k=0;k<klen;k++)
      {
        Tkvol& kvol=_kspace.getkvol(k);
        const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    VectorT3 vg=mode.getv();
	    Field& eField=mode.getfield();
	    TArray& eArray=dynamic_cast<TArray&>(eField[_cells]);
	    for(int j=0;j<neibcount;j++)
	      {
		const int f=_cellFaces(cell,j);
		if(_BCfArray[f]==3)
		  {
		    int Fgid=findFgId(f);
		    const T refl=(*(_bcMap[Fgid]))["specifiedReflection"];
		    const T oneMinusRefl=1.-refl;
		    int cell2=_faceCells(f,1);
		    VectorT3 Af=_faceArea[f];
		    
		    if(cell2==cell)
 		      {
			Af=Af*-1.;
			cell2=_faceCells(f,0);
		      }
		    
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    
		    if(VdotA<0)
		      eArray[cell2]+=mode.calce0(SumEVdotA[j])*oneMinusRefl;
		    else
		      eArray[cell2]+=mode.calce0(SumEVdotA[j])*oneMinusRefl;

		  }
	      }
	  }
      }
    */
  }

  void updateGhostCoarse(const int cell)
  {
    const int neibcount=_cellFaces.getCount(cell);
    const int klen=_kspace.getlength();
    //TArray SumEVdotA(neibcount);
    //TArray& Tl=dynamic_cast<TArray&>(_macro.temperature[_cells]);
    //T DK3=_kspace.getDK3();

    //SumEVdotA.zero();

    for(int j=0;j<neibcount;j++)
      {

	const int f=_cellFaces(cell,j);
	if(_BCfArray[f]==3)
	  {
	    int Fgid=findFgId(f);
	    const T refl=(*(_bcMap[Fgid]))["specifiedReflection"];
	    int cell2=_faceCells(f,1);
	    VectorT3 Af=_faceArea[f];
		    
	    if(cell2==cell)
	      {
		Af=Af*-1.;
		cell2=_faceCells(f,0);
	      }

	    const int cellstart=_kspace.getGlobalIndex(cell,0);
	    int cellIndex=_kspace.getGlobalIndex(cell,0);
	    int cell2Index=_kspace.getGlobalIndex(cell2,0);

	    //first loop to sum up the diffuse reflection
	    for(int k=0;k<klen;k++)
	      {
		Tkvol& kvol=_kspace.getkvol(k);
		const int numModes=kvol.getmodenum();
		//T dk3=kvol.getdk3();
		for(int m=0;m<numModes;m++)
		  {
		    Tmode& mode=kvol.getmode(m);
		    VectorT3 vg=mode.getv();
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    
		    if(VdotA>0)
		      {
			_eArray[cell2Index]=_eArray[cellIndex]*refl;/*
			Refl_pair& rpairs=mode.getReflpair(Fgid);
			const int kk=rpairs.first.second;
			Tkvol& kkvol=_kspace.getkvol(kk);
			Tmode& refMode=kkvol.getmode(m);
			const int rIndex=cellstart+refMode.getIndex()-1;
			_eArray[cell2Index]=_eArray[rIndex]*refl;*/
		      }
		    else
		      {
			Refl_pair& rpairs=mode.getReflpair(Fgid);
			const int kk=rpairs.second.second;
			Tkvol& kkvol=_kspace.getkvol(kk);
			Tmode& refMode=kkvol.getmode(m);
			const int rIndex=cellstart+refMode.getIndex()-1;
			_eArray[cell2Index]=_eArray[rIndex]*refl;
		      }
		    cellIndex++;
		    cell2Index++;
		  }
	      }
	  }
      }

    /*
    T wallTemp(Tl[cell]);
    for(int j=0;j<neibcount;j++)
      {
	const int f=_cellFaces(cell,j);
	if(_BCfArray[f]==3)
	  {
	    int cell2=_faceCells(f,1);
	    VectorT3 Af=_faceArea[f];
	    if(cell2==cell)
	      {
		Af=Af*-1.;
		cell2=_faceCells(f,0);
	      }
	    T esum=SumEVdotA[j]*DK3;
	    _kspace.calcTemp(wallTemp,esum,Af);
	    SumEVdotA[j]=wallTemp;
	    Tl[cell2]=wallTemp;
	  }
      }

    //second loop to add in the diffuse reflection
    for(int k=0;k<klen;k++)
      {
        Tkvol& kvol=_kspace.getkvol(k);
        const int numModes=kvol.getmodenum();
	for(int m=0;m<numModes;m++)
	  {
	    Tmode& mode=kvol.getmode(m);
	    VectorT3 vg=mode.getv();
	    Field& eField=mode.getfield();
	    TArray& eArray=dynamic_cast<TArray&>(eField[_cells]);
	    for(int j=0;j<neibcount;j++)
	      {
		const int f=_cellFaces(cell,j);
		if(_BCfArray[f]==3)
		  {
		    int Fgid=findFgId(f);
		    const T refl=(*(_bcMap[Fgid]))["specifiedReflection"];
		    const T oneMinusRefl=1.-refl;
		    int cell2=_faceCells(f,1);
		    VectorT3 Af=_faceArea[f];
		    
		    if(cell2==cell)
 		      {
			Af=Af*-1.;
			cell2=_faceCells(f,0);
		      }
		    
		    const T VdotA=Af[0]*vg[0]+Af[1]*vg[1]+Af[2]*vg[2];
		    
		    if(VdotA<0)
		      eArray[cell2]+=mode.calce0(SumEVdotA[j])*oneMinusRefl;
		    else
		      eArray[cell2]+=mode.calce0(SumEVdotA[j])*oneMinusRefl;

		  }
	      }
	  }
      }
    */
  }

  void correctInterface(const int cell0, TArray& Bvec)
  {
    const int klen=_kspace.gettotmodes();
    TArray Correction(klen+1);

    const int neibcount=_cellFaces.getCount(cell0);
    for(int j=0;j<neibcount;j++)
      {
	const int f=_cellFaces(cell0,j);
	if(_BCfArray[f]==4)
	  {
	    int Fgid=findFgId(f);
	    const FaceGroup& fg=_mesh.getFaceGroup(Fgid);
	    const StorageSite& faces0=fg.site;
	    const int offset=faces0.getOffset();
	    int cell1=_faceCells(f,1);
	    if(cell1==cell0)
	      cell1=_faceCells(f,0);
	    const TKConnectivity& TKC=*(_FaceToKSC.find(Fgid)->second)[f-offset];
	    TKC.multiplySelf(Bvec,Correction);
	    Bvec.zero();
	    Distribute(cell1,Correction, Bvec);
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
  const TArray& _faceAreaMag;
  const VectorT3Array& _faceArea;
  const TArray& _cellVolume;
  const VectorT3Array& _cellCoords;
  PhononMacro& _macro;
  Tkspace& _kspace;
  COMETBCMap& _bcMap;
  const IntArray& _BCArray;
  const IntArray& _BCfArray;
  T _aveResid;
  T _residChange;
  FaceToFg _fgFinder;
  COpts _options;
  const FgTKClistMap& _FaceToKSC;
  TArray& _eArray;
  TArray& _e0Array;
  TArray& _resArray;
  
};


#endif
