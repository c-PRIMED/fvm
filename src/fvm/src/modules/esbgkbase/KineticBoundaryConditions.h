// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _KINETICBOUNDARYCONDITIONS_H_
#define _KINETICBOUNDARYCONDITIONS_H_
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

#include "StressTensor.h"
template<class X,class Diag,class OffDiag>
class KineticBoundaryConditions
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
 

  KineticBoundaryConditions(const StorageSite& faces,
			    const Mesh& mesh,
			    const GeomFields& geomFields,
			    const Quadrature<X>& quadrature,
			    MacroFields& macroFields,
			    DistFunctFields<X>& dsfPtr):
			  
    _faces(faces),
    _cells(mesh.getCells()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])),
    _quadrature(quadrature), 
    _macroFields(macroFields),
    _dsfPtr(dsfPtr),
    _faceCells(mesh.getFaceCells(_faces)),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces]))
    
  {}
    KineticModelOptions<X>&   getOptions() {return _options;}  

  void applyDiffuseWallBC(int f,const VectorX3&  WallVelocity,const X& WallTemperature) const
  {
    
    const double pi=_options.pi;
    const double epsilon=_options.epsilon_ES;
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); ///changed
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
   
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
    const XArray& wts= dynamic_cast<const XArray&>(*_quadrature.dcxyzPtr);
    VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    //XArray& pressure  = dynamic_cast<XArray&>(_macroFields.pressure[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);   
    const X uwall = WallVelocity[0];
    const X vwall = WallVelocity[1];
    const X wwall = WallVelocity[2];
    //cout << "uwall " << uwall << endl;
    const X Twall = WallTemperature;

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
	if (c_dot_en -wallV_dot_en < T_Scalar(epsilon)) //incoming
	  {
	    Dmr = Dmr - fwall*wts[j]*(c_dot_en-wallV_dot_en);
	    incomFlux=incomFlux-dsf[c0]*wts[j]*(c_dot_en -wallV_dot_en);
	  }
	else
	  {
	    Nmr = Nmr + dsf[c0]*wts[j]*(c_dot_en -wallV_dot_en);
	  }	
      }
    const X nwall = Nmr/Dmr; // wall number density for initializing Maxwellian
    density[c1]=nwall;
    //if (c0==80)cout <<"incoming" << incomFlux <<" outgoing" <<Nmr <<endl;
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
	const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];

	if (c_dot_en-wallV_dot_en < T_Scalar(epsilon))
	  {
	    dsf[c1] = nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	    //dsf[c0]=dsf[c1];  // change value of fluid boundary-cell
	  }
	else
	  dsf[c1]=dsf[c0];  
      }  
    
   
   
  }

void applyDiffuseWallBC(const VectorX3& bVelocity,const X& bTemperature) const
  {
    for (int i=0; i<_faces.getCount();i++)
    applyDiffuseWallBC(i,bVelocity,bTemperature);
  }
  
  
 void applyDiffuseWallBC(const FloatValEvaluator<VectorX3>& bVelocity,const FloatValEvaluator<X>& bTemperature)const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyDiffuseWallBC(i,bVelocity[i],bTemperature[i]);

  }

  // Real wall with momentum accomodation coefficient 
 void applyRealWallBC(int f,const VectorX3&  WallVelocity,const X& WallTemperature,const X& accommodationCoefficient,const vector<int>& vecReflection) const
  {
    
    const double pi=_options.pi;
    //const double epsilon=_options.epsilon_ES;
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); 
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
   
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
	    incomFlux=incomFlux-dsf[c0]*wts[j]*(c_dot_en-wallV_dot_en);
	  }
	else
	  {
	   Nmr = Nmr + dsf[c0]*wts[j]*(c_dot_en-wallV_dot_en);
	  }	
      }
    const X nwall = Nmr/Dmr; // wall number density for initializing Maxwellian
    density[c1]=nwall;
    //if (c0==80)cout <<"incoming " << incomFlux <<" outgoing " <<Nmr << " dens " << nwall<< endl;
    //update f in boundary cell =alpha*f_wall+(1-alpha)*f_reflection
    //c.en-vwall.en=vwall.en-c_incident.en

    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
	const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];

	// needed for method 1
	/*
	const int N2 = _quadrature.getNthetaCount();
	const int N3 = _quadrature.getNphiCount();
	const X dcx =  _quadrature.get_dcx();
	const X dcy = _quadrature.get_dcy();
	const X dcz = _quadrature.get_dcz();
	*/
	if (c_dot_en-wallV_dot_en < T_Scalar(0.0))
	  { 
	    /*
	    const X cx_incident = cx[j] - 2.0*(c_dot_en-wallV_dot_en)*en[0];
	    const X cy_incident = cy[j] - 2.0*(c_dot_en-wallV_dot_en)*en[1];
	    const X cz_incident = cz[j] - 2.0*(c_dot_en-wallV_dot_en)*en[2];    	   
	    
	    //method 1	    
	    const X i_incident = int((cx_incident-cx[0])/dcx+0.5);
	    const X j_incident = int((cy_incident-cy[0])/dcy+0.5);
	    const X k_incident = int((cz_incident-cz[0])/dcz+0.5);
	    const int direction_incident = k_incident+N3*j_incident+N3*N2*i_incident;
	    
	    //method 2
	    int direction_incident=0; 
	    X Rdotprod=1e54;
	    X dotprod=0.0;
	    for (int js=0; js<numDirections; js++){
	      dotprod=pow(cx_incident-cx[js],2)+pow(cy_incident-cy[js],2)+pow(cz_incident-cz[js],2);
		if (dotprod< Rdotprod){
		  Rdotprod =dotprod;
		  direction_incident=js;}
	    }
	   
	    */
	    const int direction_incident = vecReflection[j];
	    Field& fndi = *_dsfPtr.dsf[direction_incident];
	    const XArray& dsfi = dynamic_cast<const XArray&>(fndi[_cells]);
	    
	    dsf[c1] = alpha*nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall)+m1alpha*dsfi[c0]; //write into boundary;
	    //dsf[c0]=dsf[c1]; //write into boundary cell
	  }
	else
	  dsf[c1]=dsf[c0];  
      }  
    
   
   
  }

 void applyRealWallBC(const VectorX3& bVelocity,const X& bTemperature,const X& accomCoeff,const vector<int>& vecReflection) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyRealWallBC(i,bVelocity,bTemperature,accomCoeff,vecReflection);
  }
  
 
  void applyRealWallBC(const FloatValEvaluator<VectorX3>& bVelocity,const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& accomCoeff,const vector<int>& vecReflection) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyRealWallBC(i,bVelocity[i],bTemperature[i],accomCoeff[i],vecReflection);

  }
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  /*
  void applyInsulatedWallBC(int f,const VectorX3&  WallVelocity) const 
  {

    const double pi=_options.pi;
    const double epsilon=_options.epsilon_ES;
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); ///changed
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
   
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
    const XArray& wj= dynamic_cast<const XArray&>(*_quadrature.dcxyzPtr);
    
    VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    //XArray& pressure  = dynamic_cast<XArray&>(_macroFields.pressure[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);   
    
    const X uwall=WallVelocity[0];
    const X vwall=WallVelocity[1];
    const X wwall=WallVelocity[2];
    v[c1][0]=uwall;
    v[c1][1]=vwall;
    v[c1][2]=wwall;

    X F1 = 0.0;
    X F2 = 0.0;
    X mass_out=0.0;
    X mass_in=0.0;
    X energy_out=0.0;
    X energy_in=0.0;
    X f1_nw=0.0;
    X f1_Tw=0.0;
    X f2_nw=0.0;
    X f2_Tw=0.0;
    X nwall_i=1.0;
    X Twall_i=1.0;
 


    X Twall_Nmr=0.0; 
    X Twall_Dmr=0.0; 
    for (int d_Tw < pow(10,-13) 
	   
	   for (int j=0; j<numDirections; j++)
	     {	
	       Field& fnd = *_dsfPtr.dsf[j];
	       XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	       const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	       const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	       const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	

	       if (c_dot_en - wallV_dot_en > T_Scalar(epsilon)) //outgoing
		 {
		   mass_out = mass_out + c_dot_en * dsf[c0]*wj[j];
		   energy_out = energy_out + c_dot_en *(pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2)) * dsf[c0]*wk[j];
		 }
	       else
		 {
		   mass_in = mass_in+ c_dot_en*nwall_i/(pow(pi*Twall_i,1.5))*exp(-(pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2))/Tw_i)*wj[j];

		   energy_in=energy_in + c_dot_en*(pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2))*nwall_i/(pow(pi*Twall_i,1.5))*exp(-(pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2))/Tw_i)*wj[j];

		   f1_nw = f1_nw + c_dot_en / pow(pi*Tw_i,1.5) * exp(- (pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2))/Tw_i)*wj[j];

		   f2_nw = f2_nw + c_dot_en * (pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2)) / (pow(pi*Tw_i),1.5) * exp(-((pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2))/Tw_i)*wj[j];
																    f1_Tw_a = f1_Tw_a -3/2*c_dot_en*(pow(cx[j]-uwall,2)+pow(cy[j]-vwall)+pow(cz[j]-wwall,2))*nwall_i* wj[j]/(pow(pi,1.5)*pow(Twall_i,2.5))*exp(-(pow(cx[j]-uwall,2)+pow(cy[j]-vwall,2)+pow(cz[j]-wwall,2))/Twall_i);
		   f1_Tw_b =  f1_Tw_b +
		   f2_Tw_a =
		   f2_Tw_b =
		   
		 }

  

	
	if (c_dot_en -wallV_dot_en > T_Scalar(epsilon)) 
	  {
	    Twall_Nmr = Twall_Nmr+c_dot_en*(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))*dsf[c0]*wj[j]; 
	    Twall_Dmr = Twall_Dmr+c_dot_en*dsf[c0]*wj[j];
	  }
	
      }
    const X Twall = Twall_Nmr/(2*Twall_Dmr);
  
    temperature[c1]=Twall;
    //X Nmr(0.0) ;
    //X Dmr(0.0) ;
    //X incomFlux(0.0);
    X nw_coef = 2.0*pow(pi,0.5)/(pow(Twall,0.5));
    X nw_sum=0.0;
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	//const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	//const X fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	if (c_dot_en > T_Scalar(epsilon)) 
	  {
	    nw_sum=nw_sum+(c_dot_en)*dsf[c0]*wj[j];
	      //Dmr = Dmr - fwall*wj[j]*(c_dot_en-wallV_dot_en);
	      //incomFlux=incomFlux-dsf[c0]*wj[j]*(c_dot_en -wallV_dot_en);
	  }	
      }
    const X nwall = nw_coef*nw_sum; // wall number density for initializing Maxwellian
   

 density[c1]=nwall;
    //if (c0==80)cout <<"incoming" << incomFlux <<" outgoing" <<Nmr <<endl;
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
	const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];

	if (c_dot_en-wallV_dot_en < T_Scalar(epsilon))
	  {
	   
	    dsf[c1] = (nwall/pow(pi*Twall,1.5))*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	    //dsf[c0]=dsf[c1];  // change value of fluid boundary-cell
	  }
	else
	  dsf[c1]=dsf[c0];  
      }  

  }

   void applyInsulatedWallBC(const FloatValEvaluator<VectorX3>& bVelocity) const 
//,const FloatValEvaluator<X>& bTemperature,const vector<int>& vecReflection) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyInsulatedWallBC(i, bVelocity[i]);//,bTemperature[i]);
  }

  void applyInsulatedWallBC(const VectorX3& bVelocity)const
//, const X& bTemperature) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyInsulatedWallBC(i, bVelocity);//,bTemperature[i]);
  }
 
  
   
  */

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void applySpecularWallBC(int f,const vector<int>& vecReflection) const
  {
    
    const int c0 = _faceCells(f,0); //interior
    const int c1 = _faceCells(f,1); //boundary cell
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
    
    for (int j=0; j<numDirections; j++)
      {
	// Find incident molecule directions
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X  c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	
	Field& fndw = *_dsfPtr.dsf[j];
	XArray& dsfw = dynamic_cast<XArray&>(fndw[_cells]);
	if(c_dot_en < T_Scalar(0.0)) //incoming molecules to interior
	  {
	    const int direction_incident = vecReflection[j];
	    Field& fnd = *_dsfPtr.dsf[direction_incident];
	    const XArray& dsf = dynamic_cast<const XArray&>(fnd[_cells]);
	    dsfw[c1] = dsf[c0]; //write into boundary
	  }
	else{
	  dsfw[c1]=dsfw[c0];}
      }
    
    
  }   
  
  void applySpecularWallBC(const vector<int>& vecReflection) const
  {    
    for (int i=0; i<_faces.getCount();i++){
      applySpecularWallBC(i,vecReflection);  
      // applySpecularWallBC_Cartesian(i);
    }
    
  }
  
  void applySpecularWallBC_Cartesian(int f) const
  {
    
    const int c0 = _faceCells(f,0); //interior
    const int c1 = _faceCells(f,1); //boundary cell

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
   
    const int N2 = _quadrature.getNthetaCount();
    const int N3 = _quadrature.getNphiCount();
    const X dcx =  _quadrature.get_dcx();
    const X dcy = _quadrature.get_dcy();
    const X dcz = _quadrature.get_dcz();
   
    for (int j=0; j<numDirections; j++)
      {
	// Find incident molecule directions
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X  c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	
	Field& fndw = *_dsfPtr.dsf[j];
	XArray& dsfw = dynamic_cast<XArray&>(fndw[_cells]);
	if(c_dot_en < T_Scalar(0.0)) //incoming molecules to interior
	  {
	   
	    const X cx_incident = cx[j] - 2.0*c_dot_en*en[0];
	    const X cy_incident = cy[j] - 2.0*c_dot_en*en[1];
	    const X cz_incident = cz[j] - 2.0*c_dot_en*en[2];    	   
	    
	    const X i_incident = int((cx_incident-cx[0])/dcx+0.5);
	    const X j_incident = int((cy_incident-cy[0])/dcy+0.5);
	    const X k_incident = int((cz_incident-cz[0])/dcz+0.5);
	    const int direction_incident = k_incident+N3*j_incident+N3*N2*i_incident;
	   	    
	    Field& fnd = *_dsfPtr.dsf[direction_incident];
	    const XArray& dsf = dynamic_cast<const XArray&>(fnd[_cells]);
	    
	    dsfw[c1] = dsf[c0]; //write into boundary
	    
	  }
	else{
	  dsfw[c1]=dsfw[c0];}
      }
    
        
  }   

  void applyZeroGradientBC(int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const int numDirections = _quadrature.getDirCount();
    
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	dsf[c1]=dsf[c0];
      }
  }
  void applyZeroGradientBC() const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyZeroGradientBC(i);
  } 
  
  void applyPressureInletBC(int f,const X& inletTemperature,const X& inletPressure)  const
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
    const XArray& wts = dynamic_cast<const XArray&>(*_quadrature.dcxyzPtr);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    XArray& pressure  = dynamic_cast<XArray&>(_macroFields.pressure[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);   
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
	else{dsf[c1]=dsf[c0];}//outgoing
	
      }
    
}
  

  
  
  void applyPressureInletBC(const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& bPresssure) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyPressureInletBC(i,bTemperature[i],bPresssure[i]);

  } 
  

  // velocity inlet given temperature and mass flow rate
     
 void applyInletBC(int f,const VectorX3&  InletVelocity,const X& InletTemperature,const X& InletMdot,const vector<int>& vecReflection) const
  {
    
    const double pi=_options.pi;
    //const double epsilon=_options.epsilon_ES;
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); 
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
   
    const int numDirections = _quadrature.getDirCount();
    const XArray& cx = dynamic_cast<const XArray&>(*_quadrature.cxPtr);
    const XArray& cy = dynamic_cast<const XArray&>(*_quadrature.cyPtr);
    const XArray& cz = dynamic_cast<const XArray&>(*_quadrature.czPtr);
    const XArray& wts= dynamic_cast<const XArray&>(*_quadrature.dcxyzPtr);
    VectorX3Array& v = dynamic_cast<VectorX3Array&>(_macroFields.velocity[_cells]);
    XArray& density  = dynamic_cast<XArray&>(_macroFields.density[_cells]); 
    XArray& temperature  = dynamic_cast<XArray&>(_macroFields.temperature[_cells]);   
    const X uin = InletVelocity[0];
    const X vin = InletVelocity[1];
    const X win = InletVelocity[2];
    const X Tin = InletTemperature;
    const X mdot = InletMdot;
    
    v[c1][0]=uin;
    v[c1][1]=vin;
    v[c1][2]=win;
  
    temperature[c1]=Tin;
    X Nmr(0.0) ;
    X Dmr(0.0);	
    const X rho_init=_options["rho_init"]; 
    const X T_init= _options["T_init"]; 
    const X R=8314.0/_options["molecularWeight"];
    const X u_init=pow(2.0*R*T_init,0.5);
    // find the number density corressponding to the inlet Maxwellian
    for (int j=0; j<numDirections; j++)
      {
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	const X fwall = 1.0/pow(pi*Tin,1.5)*exp(-(pow(cx[j]-uin,2.0)+pow(cy[j]-vin,2.0)+pow(cz[j]-win,2.0))/Tin);
	if (c_dot_en < T_Scalar(0.0)){Dmr = Dmr + fwall*wts[j]*c_dot_en;}
	Nmr=mdot/(rho_init*u_init);
	   
      }
    const X nin = Nmr/Dmr; // wall number density for initializing Maxwellian
    density[c1]=nin;

    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	XArray& dsf = dynamic_cast< XArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
	
	if (c_dot_en < T_Scalar(0.0))
	  { 
	    const int direction_incident = vecReflection[j];
	    Field& fndi = *_dsfPtr.dsf[direction_incident];
	    const XArray& dsfi = dynamic_cast<const XArray&>(fndi[_cells]);
	    
	    dsf[c1] = nin/pow(pi*Tin,1.5)*exp(-(pow(cx[j]-uin,2.0)+pow(cy[j]-vin,2.0)+pow(cz[j]-win,2.0))/Tin)+dsfi[c0]; //inlet Maxwellian +reflected
	  }
	else
	  dsf[c1]=dsf[c0];  
      }       
   
  }
  void applyInletBC(const VectorX3& bVelocity,const X& bTemperature,const X& Mdot,const vector<int>& vecReflection) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyInletBC(i,bVelocity,bTemperature,Mdot,vecReflection);
  }
  
 
  void applyInletBC(const FloatValEvaluator<VectorX3>& bVelocity,const FloatValEvaluator<X>& bTemperature,const FloatValEvaluator<X>& Mdot,const vector<int>& vecReflection) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyInletBC(i,bVelocity[i],bTemperature[i],Mdot[i],vecReflection);

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
  const StorageSite& _cells;
  const Array<int>& _ibType;
  const Quadrature<X>& _quadrature;
  MacroFields& _macroFields;
  const DistFunctFields<X>& _dsfPtr;
  const CRConnectivity& _faceCells;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  KineticModelOptions<X> _options;
};
  
#endif
