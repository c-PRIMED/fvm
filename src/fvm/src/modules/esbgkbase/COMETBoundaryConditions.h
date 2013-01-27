// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _COMETBOUNDARYCONDITIONS_H_
#define _COMETBOUNDARYCONDITIONS_H_
#include <stdio.h>
#include "Mesh.h"
#include <vector>
#include <cmath>
#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"
#include "GradientModel.h"
#include "FluxLimiters.h"
#include "Limiters.h"

#include "StressTensor.h"
template<class X,class Diag,class OffDiag>
class COMETBoundaryConditions
{
 public :
  
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<int> IntArray;
  
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;  
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<X> XArray;
  typedef Vector<X,3> VectorX3; //new for esbgk
  typedef Array<VectorX3> VectorX3Array;
  typedef StressTensor<X> StressTensorX6;
  typedef Array<StressTensorX6> StressTensorArray;
  typedef Gradient<T_Scalar> GradType;
  typedef Array<GradType> GradArray;
  typedef GradientModel<T_Scalar> GradModelType;
  typedef typename GradModelType::GradMatrixType GradMatrix; 

  COMETBoundaryConditions(const StorageSite& faces,
			    const Mesh& mesh,
			    const GeomFields& geomFields,
			    const Quadrature<X>& quadrature,
			    MacroFields& macroFields,
			    DistFunctFields<X>& dsfPtr):
			  
    _faces(faces),
    _allFaces(mesh.getFaces()),
    _mesh(mesh),
    _geomFields(geomFields),
    _cells(mesh.getCells()),
    _cellFaces(mesh.getCellFaces()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])),
    _quadrature(quadrature), 
    _macroFields(macroFields),
    _dsfPtr(dsfPtr),
    _faceCells(mesh.getFaceCells(_faces)),
    _allFaceCells(mesh.getAllFaceCells()),
    _cellCells(mesh.getCellCells()),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces]))
    
  {}
    COMETModelOptions<X>&   getOptions() {return _options;}  
  
    void applyPressureInletBC(int f,const X& inletTemperature,const X& inletPressure,GradMatrix& gMat)  const
    {
    
      int surface(-1);
      const double pi=_options.pi;  
      const int c0 = _faceCells(f,0);
      const int c1 = _faceCells(f,1);
    
      if (_ibType[c0] != Mesh::IBTYPE_FLUID)
	return;
    
      GradArray Grads(_quadrature.getDirCount());
      Grads.zero();

      VectorT3 Gcoeff;
      TArray limitCoeff1(_quadrature.getDirCount());
      TArray min1(_quadrature.getDirCount());
      TArray max1(_quadrature.getDirCount());

      for(int dir=0;dir<_quadrature.getDirCount();dir++)
	{
	  limitCoeff1[dir]=T_Scalar(1.e20);
	  Field& fnd = *_dsfPtr.dsf[dir];
	  TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
	  min1[dir]=dsf[c0];
	  max1[dir]=dsf[c0];
	}

    const VectorT3Array& faceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_allFaces]);
    const VectorT3Array& cellCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_cells]);   


    const int neibcount=_cellFaces.getCount(c0);
    for(int j=0;j<neibcount;j++)
      {

        const int face=_cellFaces(c0,j);

        int cell2=_allFaceCells(face,1);
        if(cell2==c0)
          cell2=_allFaceCells(face,0);

        if(cell2==c1)surface=face;
        Gcoeff=gMat.getCoeff(c0, cell2);
        for(int dir=0;dir<_quadrature.getDirCount();dir++)
          {
            Field& fnd = *_dsfPtr.dsf[dir];
            TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
            Grads[dir].accumulate(Gcoeff,dsf[cell2]-dsf[c0]);
            if(min1[dir]>dsf[cell2])min1[dir]=dsf[cell2];
            if(max1[dir]<dsf[cell2])max1[dir]=dsf[cell2];
          }
      }

    for(int j=0;j<neibcount;j++)
      {
        const int face=_cellFaces(c0,j);

        SuperbeeLimiter lf;
        for(int dir=0;dir<_quadrature.getDirCount();dir++)
          {
            Field& fnd = *_dsfPtr.dsf[dir];
            TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
            GradType& grad=Grads[dir];
            VectorT3 fVec=faceCoords[face]-cellCoords[c0];
            T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);
            computeLimitCoeff(limitCoeff1[dir], dsf[c0], SOU, min1[dir], max1[dir], lf);
          }
      }

    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr); 
    const XArray& wts = dynamic_cast<const XArray&>(*_quadrature.dcxyzPtr);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    XArray& pressure  = dynamic_cast<XArray&>(_macroFields.pressure[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);
    //XArray& heatFlux = dynamic_cast<XArray&>(_macroFields.heatFlux[_cells]);   
    VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    
   
    X uwallnew = 0.0; 
     
    const X Tin = inletTemperature;
    const X Pin = inletPressure;
    const X nin = Pin/Tin; // wall number density for initializing Maxwellian

    //store boundary values
    temperature[c1]=inletTemperature;
    pressure[c1]=inletPressure;
    density[c1]=nin;
    v[c1][0]=0.0;
    v[c1][1]=0.0;
    v[c1][2]=0.0;

    const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
    
   
    //update velocity
    for (int j=0; j<numDirections; j++)
      {
        Field& fnd = *_dsfPtr.dsf[j];
        XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
        dsf[c1] = nin/pow(pi*Tin,1.5)*exp(-(pow(cx[j]-v[c1][0],2.0)+pow(cy[j]-v[c1][1],2.0)+pow(cz[j]-v[c1][2],2.0))/Tin);

	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];

	if (abs(en[0]) == T_Scalar(1.0)) 
	  {
	    if( c_dot_en > T_Scalar(0.0))
	      {
		uwallnew=uwallnew+cx[j]*dsf[c0]*wts[j];
	      }
	    else {uwallnew=uwallnew+cx[j]*dsf[c1]*wts[j];}
	  } 
	else if (abs(en[1]) == T_Scalar(1.0)) 
	  {
	    if( c_dot_en > T_Scalar(0.0))
	      {
		uwallnew=uwallnew+cy[j]*dsf[c0]*wts[j];
	      }
	    else {uwallnew=uwallnew+cy[j]*dsf[c1]*wts[j];}
	  }  
	else if (abs(en[2]) == T_Scalar(1.0)) 
	  {
	    if( c_dot_en > T_Scalar(0.0))
	      {
		uwallnew=uwallnew+cz[j]*dsf[c0]*wts[j];
	      }
	    else {uwallnew=uwallnew+cz[j]*dsf[c1]*wts[j];}
	  }
      }
    /*
    if (abs(en[0]) == T_Scalar(1.0))     
      v[c1][0]=uwallnew/nin;
    else if (abs(en[1]) == T_Scalar(1.0)) //incoming  {vwallnew=vwallnew+cy[j]*dsf[c1]*wts[j];}
      v[c1][1]=uwallnew/nin;
    else if(abs(en[2]) == T_Scalar(1.0))     
      v[c1][2]=uwallnew/nin;//uwall=uwall+relax*(uwallnew-uwall);
    */
    

    //update f
    for (int j=0; j<numDirections; j++)
      {
        Field& fnd = *_dsfPtr.dsf[j];
        XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
        const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
        const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
        if (c_dot_en< T_Scalar(0.0))
          {
            dsf[c1] = nin/pow(pi*Tin,1.5)*exp(-(pow(cx[j]-v[c1][0],2.0)+pow(cy[j]-v[c1][1],2.0)+pow(cz[j]-v[c1][2],2.0))/Tin);
            // dsf[c0]=dsf[c1];// change value of fluid boundary-cell
          }
        else
          {
            GradType& grad=Grads[j];
            VectorT3 fVec=faceCoords[surface]-cellCoords[c0];
            VectorT3 rVec=cellCoords[c1]-cellCoords[c0];
            T_Scalar r=gMat.computeR(grad,dsf,rVec,c0,c1);
            T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);
            //dsf[c1]=dsf[c0];
            //dsf[c1]=dsf[c0]+superbee(r)*SOU;
            dsf[c1]=dsf[c0]+limitCoeff1[j]*SOU;
          }
      }
    
    }
  
    void applyPressureInletBC(const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& bPresssure) const
    {
      GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);
      for (int i=0; i<_faces.getCount();i++)
	applyPressureInletBC(i,bTemperature[i],bPresssure[i],gradMatrix);

    }
  
  //outlet Boundary Condition
  void applyPressureOutletBC(int f,const X& outletTemperature,const X& outletPressure) const
  { 
   
    const double pi=_options.pi;
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]);
    VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    XArray& pressure  = dynamic_cast<XArray&>(_macroFields.pressure[_cells]);  
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);  
   
    // X Tout = outletTemperature;
    const X Pout = outletPressure;
    X Pcell= pressure[c0];  
    X gamma =_options.SpHeatRatio;
    const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
    X nout =density[c0];
    
    if (Pcell > Pout){
      const X avel2 =gamma*Pcell/density[c0]; //velcotiy of sound square
      nout =density[c0]-(Pcell-Pout)/avel2;
      X ubulk=pow(pow(v[c0][0],2.0)+pow(v[c0][1],2.0)+pow(v[c0][2],2.0),0.5);
      X ucoeff=1.0/ubulk*(ubulk+(Pcell-Pout)/(pow(2.0*avel2,0.5)*density[c0]));
      if(abs(en[0])==T_Scalar(1.0))
	v[c1][0] = v[c0][0]*ucoeff; 
      //v[c1][0]=v[c0][0]+(pressure[c0]-Pout)/(pow(2.0*avel2,0.5)*density[c0]); //exit velocity
      else if(abs(en[1])==T_Scalar(1.0))
	v[c1][1] = v[c0][1]*ucoeff;
      else if(abs(en[2])==T_Scalar(1.0))
	v[c1][2] = v[c0][2]*ucoeff;
      
      density[c1]=nout; pressure[c1]=Pout; //update macroPr
      temperature[c1]=Pout/nout;
    }
    else{
      
      density[c1]=density[c0];
      pressure[c1]=pressure[c0];
      temperature[c1]=temperature[c0];
      
    }
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);

	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
	if (c_dot_en < T_Scalar(0.0))
	  {
	    dsf[c1] = nout/pow(pi*temperature[c1],1.5)*exp(-(pow(cx[j]-v[c1][0],2.0)+pow(cy[j]-v[c1][1],2.0)+pow(cz[j]-v[c1][2],2.0))/temperature[c1]);
	    // dsf[c0]=dsf[c1];// change value of fluid boundary-cell
	  }
	else{dsf[c1]=dsf[c0];}
	
      }
    
    
}

  void applyPressureOutletBC(const X& bTemperature, const X& bPressure) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyPressureOutletWallBC(i,bTemperature,bPressure);
  }
  
  
  void applyPressureOutletBC(const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& bPresssure) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyPressureOutletBC(i,bTemperature[i],bPresssure[i]);

  }

  // Real wall with momentum accomodation coefficient 
  void applyRealWallBC(int f,const VectorX3&  WallVelocity,const X& WallTemperature,const X& accommodationCoefficient,const vector<int>& vecReflection,GradMatrix& gMat) const
  {
    int surface(-1);    
    const double pi=_options.pi;
    //const double epsilon=_options.epsilon_ES;
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); 
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    GradArray Grads(_quadrature.getDirCount());
    Grads.zero();

    VectorT3 Gcoeff;
    TArray limitCoeff1(_quadrature.getDirCount());
    TArray min1(_quadrature.getDirCount());
    TArray max1(_quadrature.getDirCount());
    
    for(int dir=0;dir<_quadrature.getDirCount();dir++)
      {
        limitCoeff1[dir]=T_Scalar(1.e20);
        Field& fnd = *_dsfPtr.dsf[dir];
        TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
        min1[dir]=dsf[c0];
        max1[dir]=dsf[c0];
      }

    const VectorT3Array& faceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_allFaces]);
    const VectorT3Array& cellCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_cells]);

    const int neibcount=_cellFaces.getCount(c0);
    for(int j=0;j<neibcount;j++)
      {

        const int face=_cellFaces(c0,j);

        int cell2=_allFaceCells(face,1);
        if(cell2==c0)
          cell2=_allFaceCells(face,0);

        if(cell2==c1)surface=face;
        Gcoeff=gMat.getCoeff(c0, cell2);
        for(int dir=0;dir<_quadrature.getDirCount();dir++)
          {
            Field& fnd = *_dsfPtr.dsf[dir];
            TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
            Grads[dir].accumulate(Gcoeff,dsf[cell2]-dsf[c0]);
            if(min1[dir]>dsf[cell2])min1[dir]=dsf[cell2];
            if(max1[dir]<dsf[cell2])max1[dir]=dsf[cell2];
          }
      }
    
    for(int j=0;j<neibcount;j++)
      {
        const int face=_cellFaces(c0,j);

        SuperbeeLimiter lf;
        for(int dir=0;dir<_quadrature.getDirCount();dir++)
          {
            Field& fnd = *_dsfPtr.dsf[dir];
            TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
            GradType& grad=Grads[dir];
            VectorT3 fVec=faceCoords[face]-cellCoords[c0];
            T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);
            computeLimitCoeff(limitCoeff1[dir], dsf[c0], SOU, min1[dir], max1[dir], lf);
          }
      }

    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
    const XArray& wts= dynamic_cast<const XArray&>(*_quadrature.dcxyzPtr);
    VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);   
    const X uwall = WallVelocity[0];
    const X vwall = WallVelocity[1];
    const X wwall = WallVelocity[2];
    const X Twall = WallTemperature;
    const X alpha = accommodationCoefficient;
    X m1alpha=1.0-alpha;
    v[c1][0]=uwall;
    v[c1][1]=vwall;
    v[c1][2]=wwall;
  
    temperature[c1]=Twall;
    X Nmr(0.0) ;
    X Dmr(0.0) ;
    X incomFlux(0.0);
    for (int j=0; j<numDirections; j++)
      {
        Field& fnd = *_dsfPtr.dsf[j];
        XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
        const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
        const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
        const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
        const X fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
        if (c_dot_en -wallV_dot_en < T_Scalar(0.0)) //incoming
          {
            Dmr = Dmr - fwall*wts[j]*(c_dot_en-wallV_dot_en);
          }
        else
          {
            GradType& grad=Grads[j];
            VectorT3 fVec=faceCoords[surface]-cellCoords[c0];
            VectorT3 rVec=cellCoords[c1]-cellCoords[c0];
            T_Scalar r=gMat.computeR(grad,dsf,rVec,c0,c1);
            T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);
            Nmr = Nmr + (dsf[c0]+limitCoeff1[j]*SOU)*wts[j]*(c_dot_en-wallV_dot_en);
          }
      }

    const X nwall = Nmr/Dmr; // wall number density for initializing Maxwellian
    density[c1]=nwall;

    for (int j=0; j<numDirections; j++)
      {
        Field& fnd = *_dsfPtr.dsf[j];
        XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
        const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
        const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
        const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];

        if (c_dot_en-wallV_dot_en < T_Scalar(0.0))
          { 
            const int direction_incident = vecReflection[j];
            Field& fndi = *_dsfPtr.dsf[direction_incident];
            const XArray& dsfi = dynamic_cast<const XArray&>(fndi[_cells]);
            
            dsf[c1] = alpha*nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall)+m1alpha*dsfi[c0];
          }
        else
          {
            GradType& grad=Grads[j];
            VectorT3 fVec=faceCoords[surface]-cellCoords[c0];
            VectorT3 rVec=cellCoords[c1]-cellCoords[c0];
            T_Scalar r=gMat.computeR(grad,dsf,rVec,c0,c1);
            T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);
            dsf[c1]=dsf[c0]+limitCoeff1[j]*SOU;
          }
      }
  }

  void applyRealWallBC(const VectorX3& bVelocity,const X& bTemperature,const X& accomCoeff,const vector<int>& vecReflection) const
  {
    GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);
    for (int i=0; i<_faces.getCount();i++)
      applyRealWallBC(i,bVelocity,bTemperature,accomCoeff,vecReflection,gradMatrix);
  }
  
 
  void applyRealWallBC(const FloatValEvaluator<VectorX3>& bVelocity,const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& accomCoeff,const vector<int>& vecReflection) const
  {
    GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);
    for (int i=0; i<_faces.getCount();i++)
      applyRealWallBC(i,bVelocity[i],bTemperature[i],accomCoeff[i],vecReflection,gradMatrix);

  }

  void applyZeroGradientBC(int f, GradMatrix& gMat) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;

    int surface(-1);
    
    GradArray Grads(_quadrature.getDirCount());
    Grads.zero();

    VectorT3 Gcoeff;
    TArray limitCoeff1(_quadrature.getDirCount());
    TArray min1(_quadrature.getDirCount());
    TArray max1(_quadrature.getDirCount());
    
    for(int dir=0;dir<_quadrature.getDirCount();dir++)
      {
        limitCoeff1[dir]=T_Scalar(1.e20);
        Field& fnd = *_dsfPtr.dsf[dir];
        TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
        min1[dir]=dsf[c0];
        max1[dir]=dsf[c0];
      }

    const VectorT3Array& faceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_allFaces]);
    const VectorT3Array& cellCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_cells]);

    const int neibcount=_cellFaces.getCount(c0);
    for(int j=0;j<neibcount;j++)
      {

        const int face=_cellFaces(c0,j);
        
        int cell2=_allFaceCells(face,1);
        if(cell2==c0)
          cell2=_allFaceCells(face,0);

        //const int cell2=_cellCells(c0,j);
        if(cell2==c1)surface=face;
        Gcoeff=gMat.getCoeff(c0, cell2);
        for(int dir=0;dir<_quadrature.getDirCount();dir++)
          {
            Field& fnd = *_dsfPtr.dsf[dir];
            TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
            Grads[dir].accumulate(Gcoeff,dsf[cell2]-dsf[c0]);
            if(min1[dir]>dsf[cell2])min1[dir]=dsf[cell2];
            if(max1[dir]<dsf[cell2])max1[dir]=dsf[cell2];
          }
      }
        
    for(int j=0;j<neibcount;j++)
      {
        const int face=_cellFaces(c0,j);

        SuperbeeLimiter lf;
        for(int dir=0;dir<_quadrature.getDirCount();dir++)
          {
            Field& fnd = *_dsfPtr.dsf[dir];
            TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
            GradType& grad=Grads[dir];

            VectorT3 fVec=faceCoords[face]-cellCoords[c0];
            T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);

            computeLimitCoeff(limitCoeff1[dir], dsf[c0], SOU, min1[dir], max1[dir], lf);
          }
      }

    const int numDirections = _quadrature.getDirCount();
    
    for (int j=0; j<numDirections; j++)
      {
        Field& fnd = *_dsfPtr.dsf[j];
        XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
        GradType& grad=Grads[j];
        VectorT3 fVec=faceCoords[surface]-cellCoords[c0];
        VectorT3 rVec=cellCoords[c1]-cellCoords[c0];
        T_Scalar r=gMat.computeR(grad,dsf,rVec,c0,c1);
        T_Scalar SOU=(grad[0]*fVec[0]+grad[1]*fVec[1]+grad[2]*fVec[2]);
        dsf[c1]=dsf[c0]+limitCoeff1[j]*SOU;
      }
  }
  void applyZeroGradientBC() const
  {
    GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);
    for (int i=0; i<_faces.getCount();i++)
      applyZeroGradientBC(i,gradMatrix);
  } 
  
  void applyNSInterfaceBC( )const 
//const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& bPressure,const FloatValEvaluator<VectorX3>& bVelocity, const FloatValEvaluator<VectorX3>& bStress) const
   {
     
     for (int i=0; i<_faces.getCount();i++)
       applyNSInterfaceBC(i); //,bTemperature[i],bPressure[i],bVelocity[i],bStress[i]);


     
   } 


  
  void applyNSInterfaceBC(int f)const //X& inletTemperature,const X& inletPressure, const VectorX3& inletVelocity, const VectorX3&tau)  const
			  //StressTensor<X>&tau)  const
  {
    
   const double pi=_options.pi;  
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
   
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr); 

    //XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    //XArray& pressure  = dynamic_cast<XArray&>(_macroFields.pressure[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);   
    //VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    
    XArray& viscosity  = dynamic_cast<XArray&>(_macroFields.viscosity[_cells]); 
    
    
    VectorX3Array& IntVel = dynamic_cast<VectorX3Array&>(_macroFields.InterfaceVelocity[_faces]);
    XArray& IntPress  = dynamic_cast<XArray&>(_macroFields.InterfacePressure[_faces]); 
    //XArray& IntDens  = dynamic_cast<XArray&>(_macroFields.InterfaceDensity[_faces]); 
    StressTensorArray& tau  = dynamic_cast<StressTensorArray&>(_macroFields.InterfaceStress[_faces]); 
    
    const X Tin = temperature[c0];
    const X Pin = IntPress[f];
    const X nin = Pin/Tin; // wall number density for initializing Maxwellian

   
    const X u_in=IntVel[f][0];
    const X v_in=IntVel[f][1];
    const X w_in=IntVel[f][2];

    const VectorT3 en = _faceArea[f]/_faceAreaMag[f];  

    // Extra for CE expression
    const X rho_init=_options["rho_init"]; 
    const X T_init= _options["T_init"]; 
    const X R=8314.0/_options["molecularWeight"];
    const X nondim_length=_options["nonDimLt"];
    const X mu0=rho_init*R* T_init*nondim_length/pow(2*R* T_init,0.5);  
	
    //update f to Maxwellian
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	//const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	//const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];		
	//if (c_dot_en< T_Scalar(0.0))
	  {
	    //cout << "nin " <<nin <<endl;
	    dsf[c1] = nin/pow(pi*Tin,1.5)*exp(-(pow(cx[j]-u_in,2.0)+pow(cy[j]-v_in,2.0)+pow(cz[j]-w_in,2.0))/Tin)
	    *(1-viscosity[c0]/mu0/(nin*pow(Tin,2.0))*(pow(cx[j]-u_in,2.0)*tau[f][0]+pow(cy[j]-v_in,2.0)*tau[f][1]+pow(cz[j]-w_in,2.0)*tau[f][2])
	      +2.0*(cx[j]-u_in)*(cy[j]-v_in)*tau[f][3]+2.0*(cy[j]-v_in)*(cz[j]-w_in)*tau[f][4]+2.0*(cz[j]-w_in)*(cx[j]-u_in)*tau[f][5]); //update f to ChapmannEnskog
	  }
	//else{dsf[c1]=dsf[c0];}//outgoing
	
      } 
    //update f to ChapmannEnskog

   
//dsf[c1]=dsf[c1]*(1-2*cx[j]*cy[j]*slopeL*viscosity[c]/mu0) //velocity tangential only
//C2=pow(cx[j]-u_in,2.0)+pow(cy[j]-v_in,2.0)+pow(cz[j]-w_in,2.0)
//dsf[c1]=dsf[c1]*(1+viscosity[c1]/(mu0*Pr*nin*pow(Tin,2.0))*((cx[j]-u_in)*dtbdx+(cy[j]-v_in)*dtbdy+(cz[j]-w_in)*dtbdz)*(C2-2.5)) 
//temperature difference too
 }

 
  

 protected:

  const StorageSite& _faces;
  const StorageSite& _allFaces;
  const Mesh& _mesh;
  const GeomFields& _geomFields;
  const StorageSite& _cells;
  const CRConnectivity& _cellFaces;
  const Array<int>& _ibType;
  const Quadrature<X>& _quadrature;
  MacroFields& _macroFields;
  const DistFunctFields<X>& _dsfPtr;
  const CRConnectivity& _faceCells;
  const CRConnectivity& _allFaceCells;
  const CRConnectivity& _cellCells;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  COMETModelOptions<X> _options;
};
  
#endif
