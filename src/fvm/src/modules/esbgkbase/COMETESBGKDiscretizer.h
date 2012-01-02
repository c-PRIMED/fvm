#ifndef _COMETESBGKDISCRETIZER_H_
#define _COMETESBGKDISCRETIZER_H_

#include "Mesh.h"
#include "Quadrature.h"
#include "DistFunctFields.h"
#include "CometMatrix.h"
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
#include "SquareMatrixESBGK.h"

template<class T>
class COMETESBGKDiscretizer
{

 public:
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Vector<T_Scalar,3> VectorT3;  
  typedef Array<VectorT3> VectorT3Array;
  typedef Array<T_Scalar> TArray;
  typedef MatrixJML<T> TMatrix;
  typedef CometMatrix<T> TComet;
  typedef SquareMatrixESBGK<T> TSquareESBGK;
  typedef map<int,COMETBC<T>*> COMETBCMap;
  typedef Array<int> IntArray;
  typedef Array<bool> BoolArray;
  typedef Vector<int,2> VecInt2;
  typedef map<int,VecInt2> FaceToFg;
  typedef Vector<T,5> VectorT5; 
  typedef Array<VectorT5> VectorT5Array;
  typedef Vector<T,6> VectorT6; 
  typedef Array<VectorT6> VectorT6Array;
  typedef Vector<T,10> VectorT10; 
  typedef Array<VectorT10> VectorT10Array;
  typedef DistFunctFields<T> TDistFF;
  typedef Quadrature<T> TQuad;

 COMETESBGKDiscretizer(const Mesh& mesh, const GeomFields& geomfields, 
		       MacroFields& macroFields, TQuad& quadrature, TDistFF& dsfPtr, 
		       TDistFF& dsfPtr1, TDistFF& dsfPtr2, TDistFF& dsfEqPtrES, TDistFF& dsfPtrRes, TDistFF& dsfPtrFAS,
		       const T dT, const int order, const bool transient,
		       const T rho_init, const T T_init, const T MW,
		       COMETBCMap& bcMap, map<int, vector<int> > faceReflectionArrayMap,
		       const IntArray& BCArray, const IntArray& BCfArray,const IntArray& ZCArray):
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
    _macroFields(macroFields),
    _quadrature(quadrature),
    _dsfPtr(dsfPtr),
    _dsfPtr1(dsfPtr1),
    _dsfPtr2(dsfPtr2),
    _dsfEqPtrES(dsfEqPtrES),
    _dsfPtrRes(dsfPtrRes),
    _dsfPtrFAS(dsfPtrFAS),
    _dT(dT),
    _order(order),
    _transient(transient),
    _rho_init(rho_init),
    _T_init(T_init),
    _MW(MW),
    _bcMap(bcMap),
    _faceReflectionArrayMap(faceReflectionArrayMap),
    _BCArray(BCArray),
    _BCfArray(BCfArray),
    _ZCArray(ZCArray),
    _aveResid(-1.),
    _residChange(-1.),
    _fgFinder()
    {}
   
   void COMETSolve(const int sweep, const int level)
   {
    const int cellcount=_cells.getSelfCount();
    const int numDir = _quadrature.getDirCount();

    int start;

    if(sweep==1)
      start=0;
    if(sweep==-1)
      start=cellcount-1;

    for(int c=start;((c<cellcount)&&(c>-1));c+=sweep)
    {
        TArray Bvec(numDir+3);
	TArray Resid(numDir+3);

	if(_BCArray[c]==0)
	{
	    TComet AMat(numDir+3);   

	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    if(_transient)
	      COMETUnsteady(c,&AMat,Bvec);	
	
	    COMETConvection(c,AMat,Bvec,cellcount);
	    COMETCollision(c,&AMat,Bvec);
	    COMETMacro(c,&AMat,Bvec);

            if(level>0)
              addFAS(c,Bvec);

	    Resid=Bvec;
	    
	    AMat.Solve(Bvec);
	    Distribute(c,Bvec,Resid);
	    ComputeMacroparameters(c);
	    
	}
        else if(_BCArray[c]==1)
	{
            TSquareESBGK AMat(numDir+3);

            Bvec.zero();
            Resid.zero();
            AMat.zero();

            if(_transient)
              COMETUnsteady(c,&AMat,Bvec);

            COMETConvection(c,AMat,Bvec);
            COMETCollision(c,&AMat,Bvec);
            COMETMacro(c,&AMat,Bvec);

            if(level>0)
              addFAS(c,Bvec);
           
	    Resid=Bvec;

	    AMat.Solve(Bvec);
	    Distribute(c,Bvec,Resid);
	    ComputeMacroparameters(c);
	    
	}
        else
          throw CException("Unexpected value for boundary cell map.");
    }
  }

  void COMETUnsteady(const int cell, TMatrix* Amat, TArray& BVec)
  {
    const int numDir = _quadrature.getDirCount();
    const T two(2.0);
    const T onePointFive(1.5);
    const T pointFive(0.5);

    int count = 1;

    for(int direction=0;direction<numDir;direction++)
    {
	const T fbydT = _cellVolume[cell]/_dT; //pow(_nonDimLength,3);
        Field& fnd = *_dsfPtr.dsf[direction];
        const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
        Field& fN1nd = *_dsfPtr1.dsf[direction];
        const TArray& fN1 = dynamic_cast<const TArray&>(fN1nd[_cells]);
        Field& fN2nd = *_dsfPtr2.dsf[direction];
        const TArray& fN2 = dynamic_cast<const TArray&>(fN2nd[_cells]);
	if(_order>1)
	{
	    Amat->getElement(count,count) -= fbydT*(onePointFive*f[cell]- two*fN1[cell]
						    + pointFive*fN2[cell]);
	    BVec[count-1] -= fbydT*onePointFive;
	}
	else
	{
	    Amat->getElement(count,count) -= fbydT;
	    BVec[count-1] -= fbydT*(f[cell]- fN1[cell]);
	}
	count++;
    }
  }

  void COMETConvection(const int cell, TComet& Amat, TArray& BVec, const int cellcount)
  {
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
    const int numDir = _quadrature.getDirCount();

    const int neibcount=_cellFaces.getCount(cell);

    for(int j=0;j<neibcount;j++)
    {
        const int f=_cellFaces(cell,j);
        int cell2=_faceCells(f,1);
        VectorT3 Af=_faceArea[f];
        VectorT3 en = _faceArea[f]/_faceAreaMag[f];

        T flux;

        if(cell2==cell)
	{
            Af=Af*(-1.);
            en=en*(-1.);
            cell2=_faceCells(f,0);
	}
	
	int count=1;
	for(int dir=0;dir<numDir;dir++)
	{
	    Field& fnd = *_dsfPtr.dsf[dir];
	    const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
	    flux=cx[dir]*Af[0]+cy[dir]*Af[1]+cz[dir]*Af[2];
	    const T c_dot_en = cx[dir]*en[0]+cy[dir]*en[1]+cz[dir]*en[2];
	    
	    if(c_dot_en>T_Scalar(0))
	    {
		Amat.getElement(count,count)-=flux;
		BVec[count-1]-=flux*f[cell];
	    }
	    else
	    {
		if(_ZCArray[cell]==1)
		{
		    if(cell2>=cellcount)
		    {
			Amat.getElement(count,count)-=flux;
			BVec[count-1]-=flux*f[cell];
		    }
		    else
		      BVec[count-1]-=flux*f[cell2];
		  
		} 
		else
		  BVec[count-1]-=flux*f[cell2];
	    }
	    count++;
	}
    }
  }


  void COMETConvection(const int cell, TSquareESBGK& Amat, TArray& BVec)
  {
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);
    const int numDir = _quadrature.getDirCount();

    const int neibcount=_cellFaces.getCount(cell);

    const T one(1.0);

    for(int j=0;j<neibcount;j++)
    {
	const int f=_cellFaces(cell,j);
	int cell2=_faceCells(f,1);
	VectorT3 Af=_faceArea[f];
	VectorT3 en = _faceArea[f]/_faceAreaMag[f];

	T flux;

	if(cell2==cell)
	{
	    Af=Af*(-1.);
	    en=en*(-1.);
	    cell2=_faceCells(f,0);
	}	
       
	if(_BCfArray[f]==2) //If the face in question is a reflecting wall
	{
	    int Fgid=findFgId(f);
	    T uwall = (*(_bcMap[Fgid]))["specifiedXVelocity"];
	    T vwall = (*(_bcMap[Fgid]))["specifiedYVelocity"];
	    T wwall = (*(_bcMap[Fgid]))["specifiedZVelocity"];
	    T Twall = (*(_bcMap[Fgid]))["specifiedTemperature"];
	    const T wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	    map<int, vector<int> >::iterator pos = _faceReflectionArrayMap.find(Fgid);
	    const vector<int>& vecReflection=(*pos).second;
	    T alpha=(*(_bcMap[Fgid]))["accommodationCoefficient"];
	    T m1alpha = one-alpha;
	    const T pi(acos(-1.0));
	    
	    //first sweep - have to calculate wall number density
	    T Nmr(0.0);
	    T Dmr(0.0);
	    int count=1;
	    for(int dir1=0;dir1<numDir;dir1++)
	    {
		Field& fnd = *_dsfPtr.dsf[dir1];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		const T fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[dir1]-uwall,2.0)+pow(cy[dir1]-vwall,2.0)+pow(cz[dir1]-wwall,2.0))/Twall);
		flux=cx[dir1]*Af[0]+cy[dir1]*Af[1]+cz[dir1]*Af[2];
		const T c_dot_en = cx[dir1]*en[0]+cy[dir1]*en[1]+cz[dir1]*en[2];
		if((c_dot_en-wallV_dot_en)>T_Scalar(0))
		{
		    Amat(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		    Nmr = Nmr + f[cell]*wts[dir1]*(c_dot_en -wallV_dot_en);
		}
		else
		{   //have to move through all other directions
		    Dmr = Dmr - fwall*wts[dir1]*(c_dot_en-wallV_dot_en);
		}
		count++;
	    }
	    
	    const T nwall = Nmr/Dmr; // wall number density for initializing Maxwellian
	    
	    //Second sweep
	    const T zero(0.0);
	    count=1;
	    for(int dir1=0;dir1<numDir;dir1++)
	    {		
	        Field& fnd = *_dsfPtr.dsf[dir1];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		flux=cx[dir1]*Af[0]+cy[dir1]*Af[1]+cz[dir1]*Af[2];
		const T c1_dot_en = cx[dir1]*en[0]+cy[dir1]*en[1]+cz[dir1]*en[2];
		if((c1_dot_en-wallV_dot_en)<T_Scalar(0))
		{
		    const T coeff1 = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[dir1]-uwall,2.0)+pow(cy[dir1]-vwall,2.0)+pow(cz[dir1]-wwall,2.0))/Twall);
		    if(m1alpha!=zero)
		    {
			const int direction_incident = vecReflection[dir1];
			Field& fndi = *_dsfPtr.dsf[direction_incident];
			const TArray& dsfi = dynamic_cast<const TArray&>(fndi[_cells]);
			Amat(count,direction_incident+1)-=flux*m1alpha;
			BVec[count-1]-=flux*m1alpha*dsfi[cell];
		    }
		    int ccount=1;
		    for(int dir2=0;dir2<numDir;dir2++)
		    {
			Field& f2nd = *_dsfPtr.dsf[dir2];
			const TArray& f2 = dynamic_cast<const TArray&>(f2nd[_cells]);
			const T c2_dot_en = cx[dir2]*en[0]+cy[dir2]*en[1]+cz[dir2]*en[2];
			if((c2_dot_en -wallV_dot_en)>T_Scalar(0))
			{
			    Field& f2nd = *_dsfPtr.dsf[dir2];
			    const TArray& f2 = dynamic_cast<const TArray&>(f2nd[_cells]);
			    T coeff2 = wts[dir2]*(c2_dot_en-wallV_dot_en);
			    Amat(count,ccount)-=flux*(coeff1*coeff2/Dmr)*alpha;
			    BVec[count-1]-=flux*(coeff1*coeff2/Dmr)*alpha*f2[cell];
			}
			ccount++;
		    }
		}
		count++;
	    } 
	}
	else if(_BCfArray[f]==3) //If the face in question is a inlet velocity face
	{ 
	    int Fgid=findFgId(f); 
	    map<int, vector<int> >::iterator pos = _faceReflectionArrayMap.find(Fgid);
	    const vector<int>& vecReflection=(*pos).second;
	    const T pi(acos(-1.0));
	    
	    //first sweep - have to calculate Dmr
	    const T uin = (*(_bcMap[Fgid]))["specifiedXVelocity"];
	    const T vin = (*(_bcMap[Fgid]))["specifiedYVelocity"];
	    const T win = (*(_bcMap[Fgid]))["specifiedZVelocity"];
	    const T Tin = (*(_bcMap[Fgid]))["specifiedTemperature"];
	    const T mdot = (*(_bcMap[Fgid]))["specifiedMassFlowRate"];
	    T Nmr(0.0);
	    T Dmr(0.0);
	    const T R=8314.0/_MW;
	    const T u_init=pow(2.0*R*_T_init,0.5);
	    int count=1;
	    for(int dir1=0;dir1<numDir;dir1++)
	    {
		Field& fnd = *_dsfPtr.dsf[dir1];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		const T fwall = 1.0/pow(pi*Tin,1.5)*exp(-(pow(cx[dir1]-uin,2.0)+pow(cy[dir1]-vin,2.0)+pow(cz[dir1]-win,2.0))/Tin);
		flux=cx[dir1]*Af[0]+cy[dir1]*Af[1]+cz[dir1]*Af[2];
		const T c_dot_en = cx[dir1]*en[0]+cy[dir1]*en[1]+cz[dir1]*en[2];
		if(c_dot_en>T_Scalar(0))
		{
		    Amat(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		}
		else
		{   //have to move through all other directions
		    Dmr = Dmr + fwall*wts[dir1]*c_dot_en;
		}
		count++;
	    }
	    
	    Nmr=mdot/(_rho_init*u_init);
	    const T nin = Nmr/Dmr; // wall number density for initializing Maxwellian
	    
	    //Second sweep
	    count=1;
	    for(int dir1=0;dir1<numDir;dir1++)
	    {
		Field& fnd = *_dsfPtr.dsf[dir1];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		flux=cx[dir1]*Af[0]+cy[dir1]*Af[1]+cz[dir1]*Af[2];
		const T c_dot_en = cx[dir1]*en[0]+cy[dir1]*en[1]+cz[dir1]*en[2];
		if(c_dot_en<T_Scalar(0))
		{
		    const int direction_incident = vecReflection[dir1];
		    Field& fndi = *_dsfPtr.dsf[direction_incident];
		    const TArray& dsfi = dynamic_cast<const TArray&>(fndi[_cells]);
		    Amat(count,direction_incident+1)-=flux;
		    BVec[count-1]-=flux*(nin/pow(pi*Tin,1.5)*exp(-(pow(cx[dir1]-uin,2.0)+pow(cy[dir1]-vin,2.0)+pow(cz[dir1]-win,2.0))/Tin)+dsfi[cell]);
		}
		count++;
	    }
	}
        else if(_BCfArray[f]==4)  //if the face in question is zero derivative
	{
            int count=1;
            for(int dir=0;dir<numDir;dir++)
	    {
                Field& fnd = *_dsfPtr.dsf[dir];
                const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
                flux=cx[dir]*Af[0]+cy[dir]*Af[1]+cz[dir]*Af[2];
                const T c_dot_en = cx[dir]*en[0]+cy[dir]*en[1]+cz[dir]*en[2];

                if(c_dot_en>T_Scalar(0))
		{
                    Amat(count,count)-=flux;
                    BVec[count-1]-=flux*f[cell];
		}
                else
		{
		    Amat(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		}
		count++;
	    }
	}
        else if(_BCfArray[f]==5)  //if the face in question is specified pressure
	{
            int count=1;
            for(int dir=0;dir<numDir;dir++)
	    {
                Field& fnd = *_dsfPtr.dsf[dir];
                const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
                flux=cx[dir]*Af[0]+cy[dir]*Af[1]+cz[dir]*Af[2];
                const T c_dot_en = cx[dir]*en[0]+cy[dir]*en[1]+cz[dir]*en[2];

                if(c_dot_en>T_Scalar(0))
		{
                    Amat(count,count)-=flux;
                    BVec[count-1]-=flux*f[cell];
		}
                else
                  BVec[count-1]-=flux*f[cell2];
                count++;
	    }
	}
	else if(_BCfArray[f]==0)  //if the face in question is not reflecting
	{
	    int count=1;
	    for(int dir=0;dir<numDir;dir++)
	    {
	        Field& fnd = *_dsfPtr.dsf[dir];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		flux=cx[dir]*Af[0]+cy[dir]*Af[1]+cz[dir]*Af[2];
		const T c_dot_en = cx[dir]*en[0]+cy[dir]*en[1]+cz[dir]*en[2];
		
		if(c_dot_en>T_Scalar(0))
		{
		    Amat(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		}
		else
		  BVec[count-1]-=flux*f[cell2];
		count++;
	    }
	}
    }
  }

  void COMETCollision(const int cell, TMatrix* Amat, TArray& BVec)
  {
    TArray& collisionFrequency = dynamic_cast<TArray&>(_macroFields.collisionFrequency[_cells]);
    const int numDir = _quadrature.getDirCount();
    const int order=numDir;

    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

    VectorT10Array& coeffg = dynamic_cast<VectorT10Array&>(_macroFields.coeffg[_cells]);
    VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[_cells]);
    
    const T two(2.0);

    T coeff;
    int count = 1;
    
    for(int direction=0;direction<numDir;direction++)
    {
        Field& fnd = *_dsfPtr.dsf[direction];
	const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]); 
	coeff =_cellVolume[cell]*collisionFrequency[cell];

	T C1=(cx[direction]-v[cell][0]);
	T C2=(cy[direction]-v[cell][1]);
	T C3=(cz[direction]-v[cell][2]);
	T fGamma=coeffg[cell][0]*exp(-coeffg[cell][1]*pow(C1,2)+coeffg[cell][2]*C1
				     -coeffg[cell][3]*pow(C2,2)+coeffg[cell][4]*C2
				     -coeffg[cell][5]*pow(C3,2)+coeffg[cell][6]*C3
				     +coeffg[cell][7]*cx[direction]*cy[direction]
				     +coeffg[cell][8]*cy[direction]*cz[direction]
				     +coeffg[cell][9]*cz[direction]*cx[direction]);
	
	Amat->getElement(count,order+1)+=coeff*fGamma*(two*coeffg[cell][1]*C1-coeffg[cell][2]);
	Amat->getElement(count,order+2)+=coeff*fGamma*(two*coeffg[cell][3]*C2-coeffg[cell][4]);
	Amat->getElement(count,order+3)+=coeff*fGamma*(two*coeffg[cell][5]*C3-coeffg[cell][6]);
	Amat->getElement(count,count)-=coeff;
	
        BVec[count-1]+=coeff*fGamma;
	BVec[count-1]-=coeff*f[cell];
	count++;
    }
  }

  void COMETMacro(const int cell, TMatrix* Amat, TArray& BVec)
  {
    const int numDir = _quadrature.getDirCount();

    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);

    const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

    VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[_cells]);

    T density(0.);
    for(int dir=0;dir<numDir;dir++)
    {
        Field& fnd = *_dsfPtr.dsf[dir];
	const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
	density+=f[cell]*wts[dir];
    }
    
    int count = 1;

    for(int dir=0;dir<numDir;dir++)
    {
        Field& fnd = *_dsfPtr.dsf[dir];
	const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
        T C1=(cx[dir]-v[cell][0]);
        T C2=(cy[dir]-v[cell][1]);
        T C3=(cz[dir]-v[cell][2]);

        Amat->getElement(numDir+1,count)+=wts[dir]*C1/density;
	BVec[numDir]+=cx[dir]*wts[dir]*f[cell]/density;
        Amat->getElement(numDir+2,count)+=wts[dir]*C2/density;
	BVec[numDir+1]+=cy[dir]*wts[dir]*f[cell]/density;
        Amat->getElement(numDir+3,count)+=wts[dir]*C3/density;
	BVec[numDir+2]+=cz[dir]*wts[dir]*f[cell]/density;
	count++;
    }
    Amat->getElement(numDir+1,numDir+1)-=1;
    BVec[numDir]-=v[cell][0];
    Amat->getElement(numDir+2,numDir+2)-=1;
    BVec[numDir+1]-=v[cell][1];
    Amat->getElement(numDir+3,numDir+3)-=1;
    BVec[numDir+2]-=v[cell][2];
  }

  void Distribute(const int cell, TArray& BVec, TArray& Rvec)
  {
    const int numDir = _quadrature.getDirCount();
    VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[_cells]);
    VectorT3Array& vR = dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[_cells]);

    for(int direction=0;direction<numDir;direction++)
    {
        Field& fnd = *_dsfPtr.dsf[direction];
	Field& fndRes = *_dsfPtrRes.dsf[direction];
        TArray& f = dynamic_cast<TArray&>(fnd[_cells]);
	TArray& fRes = dynamic_cast<TArray&>(fndRes[_cells]);
        f[cell]-=BVec[direction];
	fRes[cell]=-Rvec[direction];
    }
    v[cell][0]-=BVec[numDir];
    v[cell][1]-=BVec[numDir+1];
    v[cell][2]-=BVec[numDir+2];
    vR[cell][0]=-Rvec[numDir];
    vR[cell][1]=-Rvec[numDir+1];
    vR[cell][2]=-Rvec[numDir+2];
  }

  void ComputeMacroparameters(const int cell)
  {
    const int numDir = _quadrature.getDirCount();

    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const TArray& wts= dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr);

    TArray& density = dynamic_cast<TArray&>(_macroFields.density[_cells]);
    VectorT3Array& v = dynamic_cast<VectorT3Array&>(_macroFields.velocity[_cells]);
    TArray& temperature = dynamic_cast<TArray&>(_macroFields.temperature[_cells]);
    TArray& pressure = dynamic_cast<TArray&>(_macroFields.pressure[_cells]);
    VectorT6Array& stress = dynamic_cast<VectorT6Array&>(_macroFields.Stress[_cells]);
 
    const T zero(0.);

    density[cell]=zero;
    temperature[cell]=zero;
    stress[cell][0]=0.0;stress[cell][1]=0.0;stress[cell][2]=0.0;
    stress[cell][3]=0.0;stress[cell][4]=0.0;stress[cell][5]=0.0;

    for(int dir=0;dir<numDir;dir++)
    {
        Field& fnd = *_dsfPtr.dsf[dir];
	const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
	density[cell] = density[cell]+wts[dir]*f[cell];
	temperature[cell]= temperature[cell]+(pow(cx[dir],2.0)+pow(cy[dir],2.0)
					      +pow(cz[dir],2.0))*f[cell]*wts[dir];
    }
	  	
    temperature[cell]=temperature[cell]-(pow(v[cell][0],2.0)
					 +pow(v[cell][1],2.0)
					 +pow(v[cell][2],2.0))*density[cell];
    temperature[cell]=temperature[cell]/(1.5*density[cell]);  
    pressure[cell]=density[cell]*temperature[cell];

    for(int dir=0;dir<numDir;dir++)
    {	  
	Field& fnd = *_dsfPtr.dsf[dir];
	const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);	  
	stress[cell][0] +=pow((cx[dir]-v[cell][0]),2.0)*f[cell]*wts[dir];
	stress[cell][1] +=pow((cy[dir]-v[cell][1]),2.0)*f[cell]*wts[dir];
	stress[cell][2] +=pow((cz[dir]-v[cell][2]),2.0)*f[cell]*wts[dir];
	stress[cell][3] +=(cx[dir]-v[cell][0])*(cy[dir]-v[cell][1])*f[cell]*wts[dir];
	stress[cell][4] +=(cy[dir]-v[cell][1])*(cz[dir]-v[cell][2])*f[cell]*wts[dir];
	stress[cell][5] +=(cz[dir]-v[cell][2])*(cx[dir]-v[cell][0])*f[cell]*wts[dir];	
    }
    //cout<<"density of cell "<<cell<<" is "<<density[cell]<<endl; 
  }  

  void findResid(const bool plusFAS)
  {
    const int cellcount=_cells.getSelfCount();
    const int numDir = _quadrature.getDirCount();

    TArray ResidSum(numDir+3);
    TArray Bsum(numDir+3);
    T traceSum=0.;
    T ResidScalar=0.;
    ResidSum.zero();
    Bsum.zero();
    
    for(int c=0;c<cellcount;c++)
    {
        TArray Bvec(numDir+3);
	TArray Resid(numDir+3);
	
	if(_BCArray[c]==0)
	{
	    TComet AMat(numDir+3);

	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    if(_transient)
	      COMETUnsteady(c,&AMat,Bvec);
	    
	    COMETConvection(c,AMat,Bvec,cellcount);
	    COMETCollision(c,&AMat,Bvec);
	    COMETMacro(c,&AMat,Bvec);
	 
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
	    TSquareESBGK AMat(numDir+3);

	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();

	    if(_transient)
	      COMETUnsteady(c,&AMat,Bvec);

	    COMETConvection(c,AMat,Bvec);
	    COMETCollision(c,&AMat,Bvec);
	    COMETMacro(c,&AMat,Bvec);

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
	else
          throw CException("Unexpected value for boundary cell map.");
    }
    for(int o=0;o<numDir+3;o++)
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

  void addFAS(const int c, TArray& bVec)
  {
    const int numDir = _quadrature.getDirCount();
    int count=0;
    for(int dir=0;dir<numDir;dir++)
    {
	Field& fndFAS = *_dsfPtrFAS.dsf[dir];
	TArray& fasArray=dynamic_cast<TArray&>(fndFAS[_cells]);
	bVec[count]-=fasArray[c];
	count+=1;
    }
    VectorT3Array& fasArray = dynamic_cast<VectorT3Array&>(_macroFields.velocityFASCorrection[_cells]);
    bVec[numDir]-=fasArray[c][0];
    bVec[numDir+1]-=fasArray[c][1];
    bVec[numDir+2]-=fasArray[c][2];
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
    const int numDir = _quadrature.getDirCount();
    const TArray& density = dynamic_cast<const TArray&>(_macroFields.density[_cells]);
    const VectorT3Array& v = dynamic_cast<const VectorT3Array&>(_macroFields.velocity[_cells]);
    const TArray& temperature = dynamic_cast<const TArray&>(_macroFields.temperature[_cells]);
    int count=0;
    for(int dir=0;dir<numDir;dir++)
    {
	Field& fnd = *_dsfPtr.dsf[dir];
        const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
	o[count]=f[c];
	count+=1;
    }
    o[numDir]=v[c][0];
    o[numDir+1]=v[c][1];
    o[numDir+2]=v[c][2];
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
  MacroFields& _macroFields;
  TQuad& _quadrature;
  TDistFF& _dsfPtr;
  TDistFF& _dsfPtr1;
  TDistFF& _dsfPtr2;
  TDistFF& _dsfEqPtrES;
  TDistFF& _dsfPtrRes;
  TDistFF& _dsfPtrFAS;
  const T _dT;
  const int _order;
  const bool _transient;
  const T _rho_init;
  const T _T_init;
  const T _MW;
  COMETBCMap& _bcMap;
  map<int, vector<int> > _faceReflectionArrayMap; 
  const IntArray& _BCArray;
  const IntArray& _BCfArray;
  const IntArray& _ZCArray;
  T _aveResid;
  T _residChange;
  FaceToFg _fgFinder;
  
};


#endif
