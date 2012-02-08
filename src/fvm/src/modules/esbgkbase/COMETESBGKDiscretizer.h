#ifndef _COMETESBGKDISCRETIZER_H_
#define _COMETESBGKDISCRETIZER_H_

#include "Mesh.h"
#include "Quadrature.h"
#include "DistFunctFields.h"
#include "CometMatrix.h"
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
#include "SquareMatrixESBGK.h"

#define SQR(x) ((x)*(x))

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
  typedef ArrowHeadMatrix<T,3> TArrow;
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
  _fgFinder(),
  _numDir(_quadrature.getDirCount()),
  _cx(dynamic_cast<const TArray&>(*_quadrature.cxPtr)),
  _cy(dynamic_cast<const TArray&>(*_quadrature.cyPtr)),
  _cz(dynamic_cast<const TArray&>(*_quadrature.czPtr)),
  _wts(dynamic_cast<const TArray&>(*_quadrature.dcxyzPtr)),
  _density(dynamic_cast<TArray&>(_macroFields.density[_cells])),
  _pressure(dynamic_cast<TArray&>(_macroFields.pressure[_cells])),
  _velocity(dynamic_cast<VectorT3Array&>(_macroFields.velocity[_cells])),
  _velocityResidual(dynamic_cast<VectorT3Array&>(_macroFields.velocityResidual[_cells])),
  _velocityFASCorrection(0),
  _temperature(dynamic_cast<TArray&>(_macroFields.temperature[_cells])),
  _stress(dynamic_cast<VectorT6Array&>(_macroFields.Stress[_cells])),
  _collisionFrequency(dynamic_cast<TArray&>(_macroFields.collisionFrequency[_cells])),
  _coeffg(dynamic_cast<VectorT10Array&>(_macroFields.coeffg[_cells]))

  {
    _fArrays = new TArray*[_numDir];
    _fN1Arrays = new TArray*[_numDir];
    _fN2Arrays = new TArray*[_numDir];
    _fEqESArrays = new TArray*[_numDir];
    _fResArrays = new TArray*[_numDir];
    _fasArrays = new TArray*[_numDir];
    
    for(int direction=0;direction<_numDir;direction++)
    {
        Field& fnd = *_dsfPtr.dsf[direction];
        Field& fN1nd = *_dsfPtr1.dsf[direction];
        Field& fN2nd = *_dsfPtr2.dsf[direction];
        Field& fndEqES = *_dsfEqPtrES.dsf[direction];
        Field& fndRes = *_dsfPtrRes.dsf[direction];
        Field& fndFAS = *_dsfPtrFAS.dsf[direction];

        _fArrays[direction] = &dynamic_cast<TArray&>(fnd[_cells]); 
	_fEqESArrays[direction] = &dynamic_cast<TArray&>(fndEqES[_cells]);
        _fResArrays[direction] = &dynamic_cast<TArray&>(fndRes[_cells]);
	if (fN1nd.hasArray(_cells))
	  _fN1Arrays[direction] = &dynamic_cast<TArray&>(fN1nd[_cells]);
        if (fN2nd.hasArray(_cells))
          _fN2Arrays[direction] = &dynamic_cast<TArray&>(fN2nd[_cells]);
        if (fndFAS.hasArray(_cells))
          _fasArrays[direction] = &dynamic_cast<TArray&>(fndFAS[_cells]); 
    }

    if (_macroFields.velocityFASCorrection.hasArray(_cells))
      _velocityFASCorrection = &dynamic_cast<VectorT3Array&>
        (_macroFields.velocityFASCorrection[_cells]);

  }
   void COMETSolve(const int sweep, const int level)
   {
    const int cellcount=_cells.getSelfCount();

    int start;

    if(sweep==1)
      start=0;
    if(sweep==-1)
      start=cellcount-1;

    TArray Bvec(_numDir+3);
    TArray Resid(_numDir+3);
    TArrow AMat(_numDir+3);

    for(int c=start;((c<cellcount)&&(c>-1));c+=sweep)
    {

	if(_BCArray[c]==0)
	{
	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    if(_transient)
	      COMETUnsteady(c,&AMat,Bvec);	
	
	    COMETConvection(c,AMat,Bvec,cellcount);
	    COMETTest(c,&AMat,Bvec);
	    //COMETCollision(c,&AMat,Bvec);
	    //COMETMacro(c,&AMat,Bvec);

            if(level>0)
              addFAS(c,Bvec);

	    Resid=Bvec;
	    
	    AMat.Solve(Bvec);
	    Distribute(c,Bvec,Resid);
	    //ComputeMacroparameters(c);	    
	}
        else if(_BCArray[c]==1)
	{
            Bvec.zero();
            Resid.zero();
            AMat.zero();

            if(_transient)
              COMETUnsteady(c,&AMat,Bvec);

            COMETConvection(c,AMat,Bvec);
	    COMETTest(c,&AMat,Bvec);
            //COMETCollision(c,&AMat,Bvec);
            //COMETMacro(c,&AMat,Bvec);

            if(level>0)
              addFAS(c,Bvec);
           
	    Resid=Bvec;

	    AMat.Solve(Bvec);
	    Distribute(c,Bvec,Resid);
	    //ComputeMacroparameters(c);	    
	}
        else
          throw CException("Unexpected value for boundary cell map.");
    }
  }

  template<class MatrixType>
  void COMETUnsteady(const int cell, MatrixType Amat, TArray& BVec)
  {
    const T two(2.0);
    const T onePointFive(1.5);
    const T pointFive(0.5);

    int count = 1;

    for(int direction=0;direction<_numDir;direction++)
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

  void COMETConvection(const int cell, TArrow& Amat, TArray& BVec, const int cellcount)
  {
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
	for(int dir=0;dir<_numDir;dir++)
	{
	    const TArray& f = *_fArrays[dir];
	    flux=_cx[dir]*Af[0]+_cy[dir]*Af[1]+_cz[dir]*Af[2];
	    const T c_dot_en = _cx[dir]*en[0]+_cy[dir]*en[1]+_cz[dir]*en[2];
	    
	    if(c_dot_en>T_Scalar(0))
	    {
		Amat.getElement(count,count)-=flux;
		BVec[count-1]-=flux*f[cell];
	    }
	    else
	      BVec[count-1]-=flux*f[cell2];
	    count++;
	}
    }
  }


  void COMETConvection(const int cell, TArrow& Amat, TArray& BVec)
  {

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
	    for(int dir1=0;dir1<_numDir;dir1++)
	    {
		const TArray& f = *_fArrays[dir1];
		const T fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(_cx[dir1]-uwall,2.0)+pow(_cy[dir1]-vwall,2.0)+pow(_cz[dir1]-wwall,2.0))/Twall);
		flux=_cx[dir1]*Af[0]+_cy[dir1]*Af[1]+_cz[dir1]*Af[2];
		const T c_dot_en = _cx[dir1]*en[0]+_cy[dir1]*en[1]+_cz[dir1]*en[2];
		if((c_dot_en-wallV_dot_en)>T_Scalar(0))
		{
		    Amat.getElement(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		    Nmr = Nmr + f[cell]*_wts[dir1]*(c_dot_en -wallV_dot_en);
		}
		else
		{   //have to move through all other directions
		    Dmr = Dmr - fwall*_wts[dir1]*(c_dot_en-wallV_dot_en);
		}
		count++;
	    }
	    
	    const T nwall = Nmr/Dmr; // wall number density for initializing Maxwellian

	    //Second sweep
	    const T zero(0.0);
	    count=1;
	    for(int dir1=0;dir1<_numDir;dir1++)
	    {		
		const TArray& f = *_fArrays[dir1];
		flux=_cx[dir1]*Af[0]+_cy[dir1]*Af[1]+_cz[dir1]*Af[2];
		const T c1_dot_en = _cx[dir1]*en[0]+_cy[dir1]*en[1]+_cz[dir1]*en[2];
		if((c1_dot_en-wallV_dot_en)<T_Scalar(0))
		{
		    const T coeff1 = 1.0/pow(pi*Twall,1.5)*exp(-(pow(_cx[dir1]-uwall,2.0)+pow(_cy[dir1]-vwall,2.0)+pow(_cz[dir1]-wwall,2.0))/Twall);
		    if(m1alpha!=zero)
		    {
			const int direction_incident = vecReflection[dir1];
			const TArray& dsfi = *_fArrays[direction_incident];
			//Amat.getElement(count,direction_incident+1)-=flux*m1alpha;
			BVec[count-1]-=flux*m1alpha*dsfi[cell];
		    }
		    BVec[count-1]-=flux*alpha*(nwall*coeff1);
		}
		count++;
	    } 

	}
        /*
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
	    for(int dir1=0;dir1<_numDir;dir1++)
	    {
		Field& fnd = *_dsfPtr.dsf[dir1];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		const T fwall = 1.0/pow(pi*Tin,1.5)*exp(-(pow(_cx[dir1]-uin,2.0)+pow(_cy[dir1]-vin,2.0)+pow(_cz[dir1]-win,2.0))/Tin);
		flux=_cx[dir1]*Af[0]+_cy[dir1]*Af[1]+_cz[dir1]*Af[2];
		const T c_dot_en = _cx[dir1]*en[0]+_cy[dir1]*en[1]+_cz[dir1]*en[2];
		if(c_dot_en>T_Scalar(0))
		{
		    Amat(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		}
		else
		{   //have to move through all other directions
		    Dmr = Dmr + fwall*_wts[dir1]*c_dot_en;
		}
		count++;
	    }
	    
	    Nmr=mdot/(_rho_init*u_init);
	    const T nin = Nmr/Dmr; // wall number density for initializing Maxwellian
	    
	    //Second sweep
	    count=1;
	    for(int dir1=0;dir1<_numDir;dir1++)
	    {
		Field& fnd = *_dsfPtr.dsf[dir1];
		const TArray& f = dynamic_cast<const TArray&>(fnd[_cells]);
		flux=_cx[dir1]*Af[0]+_cy[dir1]*Af[1]+_cz[dir1]*Af[2];
		const T c_dot_en = _cx[dir1]*en[0]+_cy[dir1]*en[1]+_cz[dir1]*en[2];
		if(c_dot_en<T_Scalar(0))
		{
		    const int direction_incident = vecReflection[dir1];
		    Field& fndi = *_dsfPtr.dsf[direction_incident];
		    const TArray& dsfi = dynamic_cast<const TArray&>(fndi[_cells]);
		    Amat(count,direction_incident+1)-=flux;
		    BVec[count-1]-=flux*(nin/pow(pi*Tin,1.5)*exp(-(pow(_cx[dir1]-uin,2.0)+pow(_cy[dir1]-vin,2.0)+pow(_cz[dir1]-win,2.0))/Tin)+dsfi[cell]);
		}
		count++;
	    }
	}
	*/
        else if(_BCfArray[f]==4)  //if the face in question is zero derivative
	{
            int count=1;
            for(int dir=0;dir<_numDir;dir++)
	    {
                const TArray& f = *_fArrays[dir];
                flux=_cx[dir]*Af[0]+_cy[dir]*Af[1]+_cz[dir]*Af[2];
                const T c_dot_en = _cx[dir]*en[0]+_cy[dir]*en[1]+_cz[dir]*en[2];

                if(c_dot_en>T_Scalar(0))
		{
                    Amat.getElement(count,count)-=flux;
                    BVec[count-1]-=flux*f[cell];
		}
                else
		{
		    Amat.getElement(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		}
		count++;
	    }
	}
        else if(_BCfArray[f]==5)  //if the face in question is specified pressure
	{
            int count=1;
            for(int dir=0;dir<_numDir;dir++)
	    {
                const TArray& f = *_fArrays[dir];
                flux=_cx[dir]*Af[0]+_cy[dir]*Af[1]+_cz[dir]*Af[2];
                const T c_dot_en = _cx[dir]*en[0]+_cy[dir]*en[1]+_cz[dir]*en[2];

                if(c_dot_en>T_Scalar(0))
		{
                    Amat.getElement(count,count)-=flux;
                    BVec[count-1]-=flux*f[cell];
		}
                else
                  BVec[count-1]-=flux*f[cell2];
                count++;
	    }
	}
        else if(_BCfArray[f]==6) //If the face in question is a symmetry wall
	{
            int Fgid=findFgId(f);
            map<int, vector<int> >::iterator pos = _faceReflectionArrayMap.find(Fgid);
            const vector<int>& vecReflection=(*pos).second;

            int count=1;
            for(int dir1=0;dir1<_numDir;dir1++)
	    {
                const TArray& f = *_fArrays[dir1];
                flux=_cx[dir1]*Af[0]+_cy[dir1]*Af[1]+_cz[dir1]*Af[2];
                const T c_dot_en = _cx[dir1]*en[0]+_cy[dir1]*en[1]+_cz[dir1]*en[2];
                if(c_dot_en>T_Scalar(0))
		{
                    Amat.getElement(count,count)-=flux;
                    BVec[count-1]-=flux*f[cell];
		}
                else
		{
                    const int direction_incident = vecReflection[dir1];
                    const TArray& dsfi = *_fArrays[direction_incident];
                    BVec[count-1]-=flux*dsfi[cell];
		}
                count++;
	    }
	}
	else if(_BCfArray[f]==0)  //if the face in question is not reflecting
	{
	    int count=1;
	    for(int dir=0;dir<_numDir;dir++)
	    {
                const TArray& f = *_fArrays[dir];
		flux=_cx[dir]*Af[0]+_cy[dir]*Af[1]+_cz[dir]*Af[2];
		const T c_dot_en = _cx[dir]*en[0]+_cy[dir]*en[1]+_cz[dir]*en[2];
		
		if(c_dot_en>T_Scalar(0))
		{
		    Amat.getElement(count,count)-=flux;
		    BVec[count-1]-=flux*f[cell];
		}
		else
		  BVec[count-1]-=flux*f[cell2];
		count++;
	    }
	}
    }
  }

  template<class MatrixType>
  void COMETTest(const int cell, MatrixType Amat, TArray& BVec)
  {
    const int order=_numDir;

    VectorT3Array& v = _velocity;
    TArray& density = _density;
    
    const T two(2.0);

    T coeff;
    int count = 1;
    
    for(int direction=0;direction<_numDir;direction++)
    {
	const TArray& f = *_fArrays[direction];
	const TArray& fEqES = *_fEqESArrays[direction];
	coeff =_cellVolume[cell]*_collisionFrequency[cell];
	
	T C1=(_cx[direction]-v[cell][0]);
	T C2=(_cy[direction]-v[cell][1]);
	T C3=(_cz[direction]-v[cell][2]);
	
	Amat->getElement(count,order+1)+=coeff*fEqES[cell]*(two*_coeffg[cell][1]*C1-_coeffg[cell][2]);
	Amat->getElement(count,order+2)+=coeff*fEqES[cell]*(two*_coeffg[cell][3]*C2-_coeffg[cell][4]);
	Amat->getElement(count,order+3)+=coeff*fEqES[cell]*(two*_coeffg[cell][5]*C3-_coeffg[cell][6]);
	Amat->getElement(count,count)-=coeff;
	
	BVec[count-1]+=coeff*(fEqES[cell]-f[cell]);
	
	Amat->getElement(_numDir+1,count)+=_wts[direction]*C1/density[cell];
	BVec[_numDir]+=_cx[direction]*_wts[direction]*f[cell]/density[cell];
	Amat->getElement(_numDir+2,count)+=_wts[direction]*C2/density[cell];
	BVec[_numDir+1]+=_cy[direction]*_wts[direction]*f[cell]/density[cell];
	Amat->getElement(_numDir+3,count)+=_wts[direction]*C3/density[cell];
	BVec[_numDir+2]+=_cz[direction]*_wts[direction]*f[cell]/density[cell];
	
	count++;
    }
    Amat->getElement(_numDir+1,_numDir+1)-=1;
    BVec[_numDir]-=v[cell][0];
    Amat->getElement(_numDir+2,_numDir+2)-=1;
    BVec[_numDir+1]-=v[cell][1];
    Amat->getElement(_numDir+3,_numDir+3)-=1;
    BVec[_numDir+2]-=v[cell][2];
  }
  
  template<class MatrixType>
  void COMETCollision(const int cell, MatrixType Amat, TArray& BVec)
  {
    const int order=_numDir;

    VectorT3Array& v = _velocity;
    
    const T two(2.0);

    T coeff;
    int count = 1;
    
    for(int direction=0;direction<_numDir;direction++)
    {
        const TArray& f = *_fArrays[direction];
	coeff =_cellVolume[cell]*_collisionFrequency[cell];

	T C1=(_cx[direction]-v[cell][0]);
	T C2=(_cy[direction]-v[cell][1]);
	T C3=(_cz[direction]-v[cell][2]);
	T fGamma=_coeffg[cell][0]*exp(-_coeffg[cell][1]*SQR(C1)+_coeffg[cell][2]*C1
                                      -_coeffg[cell][3]*SQR(C2)+_coeffg[cell][4]*C2
				     -_coeffg[cell][5]*SQR(C3)+_coeffg[cell][6]*C3
				     +_coeffg[cell][7]*_cx[direction]*_cy[direction]
				     +_coeffg[cell][8]*_cy[direction]*_cz[direction]
				     +_coeffg[cell][9]*_cz[direction]*_cx[direction]);
	
	Amat->getElement(count,order+1)+=coeff*fGamma*(two*_coeffg[cell][1]*C1-_coeffg[cell][2]);
	Amat->getElement(count,order+2)+=coeff*fGamma*(two*_coeffg[cell][3]*C2-_coeffg[cell][4]);
	Amat->getElement(count,order+3)+=coeff*fGamma*(two*_coeffg[cell][5]*C3-_coeffg[cell][6]);
	Amat->getElement(count,count)-=coeff;
	
        BVec[count-1]+=coeff*fGamma;
	BVec[count-1]-=coeff*f[cell];
	count++;
    }
  }

  template<class MatrixType>
  void COMETMacro(const int cell, MatrixType Amat, TArray& BVec)
  {

    VectorT3Array& v = _velocity;

    T density(0.);
    for(int dir=0;dir<_numDir;dir++)
    {
        const TArray& f = *_fArrays[dir];
	density+=f[cell]*_wts[dir];
    }
    
    int count = 1;

    for(int dir=0;dir<_numDir;dir++)
    {
        const TArray& f = *_fArrays[dir];
        T C1=(_cx[dir]-v[cell][0]);
        T C2=(_cy[dir]-v[cell][1]);
        T C3=(_cz[dir]-v[cell][2]);

        Amat->getElement(_numDir+1,count)+=_wts[dir]*C1/density;
	BVec[_numDir]+=_cx[dir]*_wts[dir]*f[cell]/density;
        Amat->getElement(_numDir+2,count)+=_wts[dir]*C2/density;
	BVec[_numDir+1]+=_cy[dir]*_wts[dir]*f[cell]/density;
        Amat->getElement(_numDir+3,count)+=_wts[dir]*C3/density;
	BVec[_numDir+2]+=_cz[dir]*_wts[dir]*f[cell]/density;
	count++;
    }
    Amat->getElement(_numDir+1,_numDir+1)-=1;
    BVec[_numDir]-=v[cell][0];
    Amat->getElement(_numDir+2,_numDir+2)-=1;
    BVec[_numDir+1]-=v[cell][1];
    Amat->getElement(_numDir+3,_numDir+3)-=1;
    BVec[_numDir+2]-=v[cell][2];
  }

  void Distribute(const int cell, TArray& BVec, TArray& Rvec)
  {
    VectorT3Array& v = _velocity;
    VectorT3Array& vR = _velocityResidual;

    for(int direction=0;direction<_numDir;direction++)
    {
	TArray& f = *_fArrays[direction];
	TArray& fRes = *_fResArrays[direction];
        f[cell]-=BVec[direction];
	fRes[cell]=-Rvec[direction];
    }
    v[cell][0]-=BVec[_numDir];
    v[cell][1]-=BVec[_numDir+1];
    v[cell][2]-=BVec[_numDir+2];
    vR[cell][0]=-Rvec[_numDir];
    vR[cell][1]=-Rvec[_numDir+1];
    vR[cell][2]=-Rvec[_numDir+2];
  }

  void ComputeMacroparameters(const int cell)
  {

    VectorT3Array& v = _velocity;
 
    const T zero(0.);

    _density[cell]=zero;
    _temperature[cell]=zero;
    _stress[cell][0]=0.0;_stress[cell][1]=0.0;_stress[cell][2]=0.0;
    _stress[cell][3]=0.0;_stress[cell][4]=0.0;_stress[cell][5]=0.0;

    for(int dir=0;dir<_numDir;dir++)
    {
	const TArray& f = *_fArrays[dir];
	_density[cell] = _density[cell]+_wts[dir]*f[cell];
	_temperature[cell]= _temperature[cell]+(SQR(_cx[dir])+SQR(_cy[dir])
					      +SQR(_cz[dir]))*f[cell]*_wts[dir];
    }
	  	
    _temperature[cell]=_temperature[cell]-(SQR(v[cell][0])
                                           +SQR(v[cell][1])
                                           +SQR(v[cell][2]))*_density[cell];
    _temperature[cell]=_temperature[cell]/(1.5*_density[cell]);  
    _pressure[cell]=_density[cell]*_temperature[cell];

    for(int dir=0;dir<_numDir;dir++)
    {	  
	const TArray& f = *_fArrays[dir];
	_stress[cell][0] +=SQR((_cx[dir]-v[cell][0]))*f[cell]*_wts[dir];
	_stress[cell][1] +=SQR((_cy[dir]-v[cell][1]))*f[cell]*_wts[dir];
	_stress[cell][2] +=SQR((_cz[dir]-v[cell][2]))*f[cell]*_wts[dir];
	_stress[cell][3] +=(_cx[dir]-v[cell][0])*(_cy[dir]-v[cell][1])*f[cell]*_wts[dir];
	_stress[cell][4] +=(_cy[dir]-v[cell][1])*(_cz[dir]-v[cell][2])*f[cell]*_wts[dir];
	_stress[cell][5] +=(_cz[dir]-v[cell][2])*(_cx[dir]-v[cell][0])*f[cell]*_wts[dir];	
    } 
  }  

  void findResid(const bool plusFAS)
  {
    const int cellcount=_cells.getSelfCount();

    TArray ResidSum(_numDir+3);
    TArray Bsum(_numDir+3);
    T traceSum=0.;
    T ResidScalar=0.;
    ResidSum.zero();
    Bsum.zero();

    TArray Bvec(_numDir+3);
    TArray Resid(_numDir+3);
    TArrow AMat(_numDir+3);
    
    for(int c=0;c<cellcount;c++)
    {
	
	if(_BCArray[c]==0)
	{
	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();
	    
	    if(_transient)
	      COMETUnsteady(c,&AMat,Bvec);
	    
	    COMETConvection(c,AMat,Bvec,cellcount);
	    COMETTest(c,&AMat,Bvec);
	    //COMETCollision(c,&AMat,Bvec);
	    //COMETMacro(c,&AMat,Bvec);
	 
            if(plusFAS)
              addFAS(c,Bvec);
   
	    traceSum+=AMat.getTraceAbs();
	    Resid+=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);

	    makeValueArray(c,Bvec);

            AMat.multiply(Resid,Bvec);
            Resid=Bvec; 

	    ArrayAbs(Bvec);
	    ArrayAbs(Resid);
	    Bsum+=Bvec;
	    ResidSum+=Resid;
	}
	else if(_BCArray[c]==1) //reflecting boundary
        {
	    Bvec.zero();
	    Resid.zero();
	    AMat.zero();

	    if(_transient)
	      COMETUnsteady(c,&AMat,Bvec);

	    COMETConvection(c,AMat,Bvec);
            COMETTest(c,&AMat,Bvec);
	    //COMETCollision(c,&AMat,Bvec);
	    //COMETMacro(c,&AMat,Bvec);

            if(plusFAS)
              addFAS(c,Bvec);

	    traceSum+=AMat.getTraceAbs();
	    Resid+=Bvec;
	    Bvec.zero();
	    Distribute(c,Bvec,Resid);

	    makeValueArray(c,Bvec);

            AMat.multiply(Resid,Bvec);
            Resid=Bvec; 

	    ArrayAbs(Bvec);
	    ArrayAbs(Resid);
	    Bsum+=Bvec;
	    ResidSum+=Resid;
	}
	else
          throw CException("Unexpected value for boundary cell map.");
    }
    for(int o=0;o<_numDir+3;o++)
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
    int count=0;
    for(int dir=0;dir<_numDir;dir++)
    {
	TArray& fasArray=*_fasArrays[dir];
	bVec[count]-=fasArray[c];
	count+=1;
    }
    VectorT3Array& fasArray = *_velocityFASCorrection;
    bVec[_numDir]-=fasArray[c][0];
    bVec[_numDir+1]-=fasArray[c][1];
    bVec[_numDir+2]-=fasArray[c][2];
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
    for(int dir=0;dir<_numDir;dir++)
    {
	const TArray& f = *_fArrays[dir];
	o[dir]=f[c];
    }
    const VectorT3Array& v = _velocity;
    o[_numDir]=v[c][0];
    o[_numDir+1]=v[c][1];
    o[_numDir+2]=v[c][2];
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
  const int _numDir;
  const TArray& _cx;
  const TArray& _cy;
  const TArray& _cz;
  const TArray& _wts;
  TArray& _density;
  TArray& _pressure;
  VectorT3Array& _velocity;
  VectorT3Array& _velocityResidual;
  VectorT3Array* _velocityFASCorrection;
  TArray& _temperature;
  VectorT6Array& _stress;
  TArray& _collisionFrequency;
  VectorT10Array& _coeffg;

  TArray** _fArrays;
  TArray** _fN1Arrays;
  TArray** _fN2Arrays;
  TArray** _fEqESArrays;
  TArray** _fResArrays;
  TArray** _fasArrays;
  
};


#endif
