// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BATTERYPCLINEARIZEINTERFACE_BV_H_
#define _BATTERYPCLINEARIZEINTERFACE_BV_H_

#include "Mesh.h"
#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "StorageSite.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "DiagonalMatrix.h"




template<class X, class Diag, class OffDiag, class otherMeshDiag>
  class BatteryPCLinearizeInterface_BV
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef Array<T_Scalar> TArray;
  typedef Array<X> XArray;
  //typedef Vector<T_Scalar,3> VectorT3;
  //typedef Array<VectorT3> VectorT3Array;
  
  //typedef DiagonalTensor<T_Scalar,2> DiagTensorT3;
  typedef CRMatrix<otherMeshDiag,otherMeshDiag,X> CCMatrix_DiagTensors;
  typedef typename CCMatrix_DiagTensors::DiagArray DiagArray_DiagTensors;

 BatteryPCLinearizeInterface_BV(const GeomFields& geomFields,
			   const T_Scalar RRConstant,
			   const T_Scalar interfaceUnderRelax,
			   const bool Anode,
			   const bool Cathode,
			   const bool bInterfaceHeatSource,
			   Field& varField):
    _geomFields(geomFields),
    _varField(varField),
    _RRConstant(RRConstant),
    _interfaceUnderRelax(interfaceUnderRelax),
    _Anode(Anode),
    _Cathode(Cathode),
    _bInterfaceHeatSource(bInterfaceHeatSource)
    {}


  void discretize(const Mesh& mesh, const Mesh& parentMesh, 
		  const Mesh& otherMesh, MultiFieldMatrix& mfmatrix, 
		  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getParentFaceGroupSite();

    // shell mesh info
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex)); 
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    const CRConnectivity& cellCells = mesh.getCellCells();
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
    DiagArray& diag = matrix.getDiag();

    // In the following, parent is assumed to be the electrolyte, and 
    // the other mesh is assumed to be electrode when implimenting
    // the Butler-Volmer equations

    // parent mesh info
    const CRConnectivity& parentFaceCells = parentMesh.getFaceCells(faces);
    const StorageSite& parentCells = parentMesh.getCells();
    const MultiField::ArrayIndex cVarIndexParent(&_varField,&parentCells);
    XArray& rParentCell = dynamic_cast<XArray&>(rField[cVarIndexParent]);
    CCMatrix_DiagTensors& parentmatrix = dynamic_cast<CCMatrix_DiagTensors&>(mfmatrix.getMatrix(cVarIndexParent,cVarIndexParent)); 
    DiagArray_DiagTensors& parentdiag = parentmatrix.getDiag();
    const TArray& faceAreaMag =
      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);

    // other mesh info
    const StorageSite& otherFaces = mesh.getOtherFaceGroupSite();
    const CRConnectivity& otherFaceCells = otherMesh.getFaceCells(otherFaces);
    const StorageSite& otherCells = otherMesh.getCells();
    const MultiField::ArrayIndex cVarIndexOther(&_varField,&otherCells);
    XArray& rOtherCell = dynamic_cast<XArray&>(rField[cVarIndexOther]);
    CCMatrix_DiagTensors& othermatrix = dynamic_cast<CCMatrix_DiagTensors&>(mfmatrix.getMatrix(cVarIndexOther,cVarIndexOther)); 
    DiagArray_DiagTensors& otherdiag = othermatrix.getDiag();

    // set constants for entire shell
    const T_Scalar F = 96485.0; //  C/mol
    const T_Scalar k = _RRConstant;
    T_Scalar csMax = 26000.0;
    if (_Anode){
      csMax = 26390.0;}
    if (_Cathode){
      csMax = 22860.0;}
    const T_Scalar alpha_a = 0.5;
    const T_Scalar alpha_c = 0.5;
    const T_Scalar R = 8.314; //  J/mol/K
    const T_Scalar dU_dT = -0.0011; // V/K
    const T_Scalar Peltier = 0.0;

    for (int f=0; f<faces.getCount(); f++)
      {
	int c0p = parentFaceCells(f,0);
	int c1p = parentFaceCells(f,1);
	if (c1p < parentCells.getSelfCount())
	  { 
	    // c0 is ghost cell and c1 is boundry cell, so swap cell numbers
	    // so that c1p refers to the ghost cell in the following
	    int temp = c0p;
	    c0p = c1p;
	    c1p = temp;
	  }

	int c0o = otherFaceCells(f,0);
	int c1o = otherFaceCells(f,1);
	if (c1o < otherCells.getSelfCount())
	  { 
	    // c0 is ghost cell and c1 is boundry cell, so swap cell numbers
	    // so that c1o refers to the ghost cell in the following
	    int temp = c0o;
	    c0o = c1o;
	    c1o = temp;
	  }

	// parent and other are full of diagonal tensors, while shell is full of square tensors
	// need to convert incoming information from parent and other
	// to a square tensor format
	OffDiag dRC0dXC3 = NumTypeTraits<OffDiag>::getUnity();
	Diag dRC0dXC0 = NumTypeTraits<Diag>::getUnity();
	OffDiag dRC0dXC2 = NumTypeTraits<OffDiag>::getUnity();
	Diag dRC0dXC1 = NumTypeTraits<Diag>::getUnity();

	// get parent and other flux
	const X parentFlux = rParentCell[c1p]; // inward shell flux on the left
	const X otherFlux = rOtherCell[c1o]; // inward shell flux on the right
	
	//get parent coeffs
	const otherMeshDiag dRC0dXC3_DiagTens = parentmatrix.getCoeff(c1p,  c0p);
	const otherMeshDiag dRC0dXC0_DiagTens = parentdiag[c1p];

	//get other coeffs
	const otherMeshDiag dRC0dXC2_DiagTens = othermatrix.getCoeff(c1o,  c0o);
	const otherMeshDiag dRC0dXC1_DiagTens = otherdiag[c1o];

	dRC0dXC3(0,0) = dRC0dXC3_DiagTens[0];
	dRC0dXC3(1,1) = dRC0dXC3_DiagTens[1];
	dRC0dXC2(0,0) = dRC0dXC2_DiagTens[0];
	dRC0dXC2(1,1) = dRC0dXC2_DiagTens[1];
	dRC0dXC1(0,0) = dRC0dXC1_DiagTens[0];
	dRC0dXC1(1,1) = dRC0dXC1_DiagTens[1];
	dRC0dXC0(0,0) = dRC0dXC0_DiagTens[0];
	dRC0dXC0(1,1) = dRC0dXC0_DiagTens[1];

	if (_bInterfaceHeatSource)
	  {
	    dRC0dXC3(2,2) = dRC0dXC3_DiagTens[2];
	    dRC0dXC2(2,2) = dRC0dXC2_DiagTens[2];
	    dRC0dXC1(2,2) = dRC0dXC1_DiagTens[2];
	    dRC0dXC0(2,2) = dRC0dXC0_DiagTens[2];
	  }

	const int c0 = f;
	const int c1 = cellCells(f,0);
	const int c2 = cellCells(f,1);
	const int c3 = cellCells(f,2);
	const T_Scalar Area = faceAreaMag[f];
	bool turnOffBV = false;

	T_Scalar cs_star = (xCell[c1])[1];
	T_Scalar ce_star = (xCell[c0])[1];
	const T_Scalar phis_star = (xCell[c1])[0];
	const T_Scalar phie_star = (xCell[c0])[0];

	// avoid nans
	if (cs_star < 0.0){ cout << "ERROR: Cs < 0, Cs=" << cs_star << endl; cs_star = 0.0;}
	if (ce_star < 0.0){ cout << "ERROR: Ce < 0, Ce=" << ce_star << endl; ce_star = 0.0;}
	//if (cs_star > csMax){ cout << "ERROR: Cs > CsMax, Cs=" << cs_star << endl; cs_star = csMax;}
	if (cs_star > csMax){ cout << "ERROR: Cs > CsMax, Cs=" << cs_star << endl; turnOffBV = true;}


	T_Scalar SOC =  cs_star/csMax;
	if (turnOffBV)
	  SOC = 1.0;

	T_Scalar U_ref = 0.1; // V
	if (_Anode){
	  U_ref = -0.16 + 1.32*exp(-3.0*SOC)+10.0*exp(-2000.0*SOC);
	  if (U_ref < 0.0)
	    {U_ref = 0.0; cout << "U_ref < 0" << endl;}
	  if (U_ref > 1.2)
	    {U_ref = 1.2; cout << "U_ref > 1.2" << endl;}
	    }
	if (_Cathode){
	  U_ref = 4.19829 + 0.0565661*tanh(-14.5546*SOC + 8.60942) - 0.0275479*(1.0/pow((0.998432-SOC),0.492465) - 1.90111) - 0.157123*exp(-0.04738*pow(SOC,8.0)) + 0.810239*exp(-40.0*(SOC-0.133875));
	  if (U_ref < 0.0)
	    {U_ref = 0.0; cout << "U_ref < 0" << endl;}
	  if (U_ref > 5.0)
	    {U_ref = 5.0; cout << "U_ref > 5.0" << endl;}
	    }

	// add small temperature dependence of U?


	T_Scalar Temp = 300.0;
	if (_bInterfaceHeatSource)
	  {
	    Temp = (xCell[c0])[2]; //  K , c0 and c1 temps are equal at convergence
	  }
	const T_Scalar C_a = alpha_a*F/R/Temp;
	const T_Scalar C_c = alpha_c*F/R/Temp;

	const T_Scalar U = U_ref - (Temp - 298.0)*dU_dT;
	const T_Scalar eta_star = phis_star - phie_star - U;

	const T_Scalar C_0 = exp(C_a*eta_star)-exp(-1.0*C_c*eta_star);

	T_Scalar i0_star = k*F*Area*pow(ce_star,alpha_c)*pow((csMax-cs_star),alpha_a)*pow(cs_star,alpha_c);
	if (turnOffBV)
	    i0_star = 0.0;

	const T_Scalar i_star = C_0*i0_star;
	
	// CURRENT SHOULD NOT BE ZERO
	//if (i_star < 1.0e-15){cout << "WARNING: current = 0" << endl;}

	// calculate dC_0/dCS
	T_Scalar dC_0dCS = 0.0;
       	const T_Scalar dC0dEta = (C_a*exp(C_a*eta_star) + C_c*exp(-1*C_c*eta_star));
	if (_Anode)
	  {
	    dC_0dCS = dC0dEta*(-1.0)*(-20000.0*exp(-2000.0*SOC) - 3.96*exp(-3.0*SOC))*(1.0/csMax);
	  }
	if (_Cathode)
	  {	 
	    dC_0dCS = dC0dEta*(-1.0)*(-0.0135664/pow((0.998432-SOC),1.49247) - 0.823297/pow(cosh(8.60942-14.5546*SOC),2.0) + 0.0595559*exp(-0.04738*pow(SOC,8.0))*pow(SOC,7.0) - 6859.94*exp(-40.0*SOC))*(1.0/csMax);
	  }

	const T_Scalar dIdCS_star = i_star*(alpha_c/cs_star - alpha_a/(csMax-cs_star)+ dC_0dCS/C_0);
	const T_Scalar dIdCE_star = i_star*alpha_c/ce_star;

	const T_Scalar dIdPhiS_star = i0_star*(C_a*exp(C_a*eta_star)+C_c*exp(-1*C_c*eta_star));
	const T_Scalar dIdPhiE_star = -1*i0_star*(C_a*exp(C_a*eta_star)+C_c*exp(-1*C_c*eta_star));

	const T_Scalar dIdTe_star = -0.5*i0_star*F*eta_star/R/Temp/Temp*(alpha_a*exp(C_a*eta_star)+alpha_c*exp(-1*C_c*eta_star));
       	const T_Scalar dIdTs_star = -0.5*i0_star*F*eta_star/R/Temp/Temp*(alpha_a*exp(C_a*eta_star)+alpha_c*exp(-1*C_c*eta_star));


	//now put flux information from meshes into shell cells

	// left(parent) shell cell - 3 neighbors
	// flux balance for all equations
	OffDiag& offdiagC0_C1 = matrix.getCoeff(c0,  c1);
	OffDiag& offdiagC0_C2 = matrix.getCoeff(c0,  c2);
	OffDiag& offdiagC0_C3 = matrix.getCoeff(c0,  c3);

	rCell[c0] = otherFlux + parentFlux;
	offdiagC0_C1 = dRC0dXC1;
	offdiagC0_C3 = dRC0dXC3;
	offdiagC0_C2 = dRC0dXC2;
	diag[c0] = dRC0dXC0;

	//Fix for species equations so that both shell cells 
	//residuals and coeffs are on same order of magnitude	
	
	(rCell[c0])[1] *= F;
	(offdiagC0_C1)(1,1) *= F;
	(offdiagC0_C3)(1,1) *= F;
	(offdiagC0_C2)(1,1) *= F;
	(diag[c0])(1,1) *= F;

	//include interface heating if thermal model turned on
	//if (_bInterfaceHeatSource)
	if (_bInterfaceHeatSource)
	  {
	    (rCell[c0])[2] += (eta_star + Peltier)*i_star; // in Watts
	    (diag[c0])(2,2) += (eta_star + Peltier)*dIdTe_star;
	    (offdiagC0_C1)(2,2) += (eta_star + Peltier)*dIdTs_star;	    
	  	
	    // Point-coupled inclusions(off diagonal terms in square tensors)
	    // Cell c0	
	    (diag[c0])(2,0) = (eta_star + Peltier)*dIdPhiE_star - i_star; 
	    (diag[c0])(2,1) = (eta_star + Peltier)*dIdCE_star; 
	    (offdiagC0_C1)(2,0) = i_star + (eta_star + Peltier)*dIdPhiS_star;
	    (offdiagC0_C1)(2,1) = i_star*dC_0dCS/dC0dEta + (eta_star + Peltier)*dIdCS_star;
	  }


	// right(other) shell cell - 2 neighbors
	// jump condition
	OffDiag& offdiagC1_C0 = matrix.getCoeff(c1,  c0);
	OffDiag& offdiagC1_C2 = matrix.getCoeff(c1,  c2);

	////////////////////////
	//       SPECIES      //
	////////////////////////

	(rCell[c1])[1] = F*otherFlux[1]  - i_star;
	offdiagC1_C0(1,1) = -1*dIdCE_star;
	offdiagC1_C2(1,1) = F*dRC0dXC2(1,1);

	//make sure diag is < 0
	if (dIdCS_star > 0.0)
	  {
	    (diag[c1])(1,1) = F*dRC0dXC1(1,1) - dIdCS_star;
	  }
	else
	  { 
	    (diag[c1])(1,1) = F*dRC0dXC1(1,1);
	  }

	////////////////////////
	//      POTENTIAL     //
	////////////////////////

	(rCell[c1])[0] = otherFlux[0] - i_star;
	offdiagC1_C0(0,0) = -1*dIdPhiE_star;
	offdiagC1_C2(0,0) = dRC0dXC2(0,0);
	(diag[c1])(0,0) = dRC0dXC1(0,0) - dIdPhiS_star;

	//make sure diag is < 0
	if ((diag[c1])(0,0) > 0.0)
	  {
	    cout << "Warning: Potential Diag > 0" << endl;
	  }

	////////////////////////
	//      THERMAL       //
	////////////////////////

	if (_bInterfaceHeatSource)
	  {
	    T_Scalar Factor = 1.0e-12; 
	    (rCell[c1])[2] = Factor*((xCell[c0])[2] - (xCell[c1])[2]);
	    offdiagC1_C0(2,2) = Factor;
	    offdiagC1_C2(2,2) = 0.0;
	    (diag[c1])(2,2) = -1.0*Factor;
	  }


	// Point-coupled inclusions(off diagonal terms in square tensors)
	// Cell c1
	(diag[c1])(0,1) = -1*dIdCS_star;
	(diag[c1])(1,0) = -1*dIdPhiS_star;
	offdiagC1_C0(0,1) = -1*dIdCE_star;
	offdiagC1_C0(1,0) = -1*dIdPhiE_star;
	
	if (_bInterfaceHeatSource)
	  {	   
	    (diag[c1])(0,2) = -1*dIdTs_star;	    
	    (diag[c1])(1,2) = -1*dIdTs_star;
	    offdiagC1_C0(0,2) = -1*dIdTe_star;
	    offdiagC1_C0(1,2) = -1*dIdTe_star;
	  }
	
	// set other coeffs to zero for right shell cell
	OffDiag& offdiagC1_C3 = matrix.getCoeff(c1,  c3);
	offdiagC1_C3 = NumTypeTraits<OffDiag>::getZero();

	// underrelax diagonal to help convergence
	diag[c1] = 1.0/_interfaceUnderRelax*diag[c1];

	if (turnOffBV)
	  {
	    cout << "r[c1]: " << (rCell[c1])[0] << (rCell[c1])[1] << endl;
	    cout << "r[c0]: " << (rCell[c0])[0] << (rCell[c1])[1] << endl;
	  }

	if (c0 == 0)
	  {

	    //cout << "OtherHeatFlux: " << otherFlux[2] << endl;
	    //cout << "ParentHeatFlux: " << parentFlux[2] << endl;
	    //cout << "SourceHeatFlux: " << (eta_star + Peltier)*i_star << endl;
	    if (_Cathode)
	      {
		//cout << "UrefAnode: " << U_ref << endl;
		cout << "Cs0: " << cs_star << endl;
		//cout << "Diag: " << diag[c1] << endl;
		//cout << "OffDiag: " << offdiagC1_C2 << endl;

		/*
		cout << " " << endl;
		cout << "OffDiag_C1C0: " << offdiagC1_C0(0,0) << " " << offdiagC1_C0(0,1) << " " << offdiagC1_C0(0,2) << endl;
		cout << "OffDiag_C1C0: " << offdiagC1_C0(1,0) << " " << offdiagC1_C0(1,1) << " " << offdiagC1_C0(1,2) << endl;
		cout << "OffDiag_C1C0: " << offdiagC1_C0(2,0) << " " << offdiagC1_C0(2,1) << " " << offdiagC1_C0(2,2) << endl;
		cout << " " << endl;
		cout << "OffDiag_C0C1: " << offdiagC0_C1(0,0) << " " << offdiagC0_C1(0,1) << " " << offdiagC0_C1(0,2) << endl;
		cout << "OffDiag_C0C1: " << offdiagC0_C1(1,0) << " " << offdiagC0_C1(1,1) << " " << offdiagC0_C1(1,2) << endl;
		cout << "OffDiag_C0C1: " << offdiagC0_C1(2,0) << " " << offdiagC0_C1(2,1) << " " << offdiagC0_C1(2,2) << endl;
		cout << " " << endl;
		cout << "Diag_C0: " << (diag[c0])(0,0) << " " << (diag[c0])(0,1) << " " << (diag[c0])(0,2) << endl;
		cout << "Diag_C0: " << (diag[c0])(1,0) << " " << (diag[c0])(1,1) << " " << (diag[c0])(1,2) << endl;
		cout << "Diag_C0: " << (diag[c0])(2,0) << " " << (diag[c0])(2,1) << " " << (diag[c0])(2,2) << endl;
		cout << " " << endl;
		cout << "Diag_C1: " << (diag[c1])(0,0) << " " << (diag[c1])(0,1) << " " << (diag[c1])(0,2) << endl;
		cout << "Diag_C1: " << (diag[c1])(1,0) << " " << (diag[c1])(1,1) << " " << (diag[c1])(1,2) << endl;
		cout << "Diag_C1: " << (diag[c1])(2,0) << " " << (diag[c1])(2,1) << " " << (diag[c1])(2,2) << endl;
		cout << " " << endl;
		cout << "Residual_C0: " << (rCell[c0])[0] << endl;
		cout << "Residual_C0: " << (rCell[c0])[1] << endl;
		cout << "Residual_C0: " << (rCell[c0])[2] << endl;
		cout << " " << endl;
		cout << "Residual_C1: " << (rCell[c1])[0] << endl;
		cout << "Residual_C1: " << (rCell[c1])[1] << endl;
		cout << "Residual_C1: " << (rCell[c1])[2] << endl;
		cout << " " << endl;
		*/
	      }
	    if (_Anode)
	      {
		//cout << "UrefAnode: " << U_ref << endl;
		//cout << "Cs0: " << cs_star << endl;
		//cout << "Diag: " << diag[c1] << endl;
		//cout << "OffDiag: " << offdiagC1_C2 << endl;

		/*
		cout << " " << endl;
		cout << "OffDiag_C1C0: " << offdiagC1_C0(0,0) << " " << offdiagC1_C0(0,1) << " " << offdiagC1_C0(0,2) << endl;
		cout << "OffDiag_C1C0: " << offdiagC1_C0(1,0) << " " << offdiagC1_C0(1,1) << " " << offdiagC1_C0(1,2) << endl;
		cout << "OffDiag_C1C0: " << offdiagC1_C0(2,0) << " " << offdiagC1_C0(2,1) << " " << offdiagC1_C0(2,2) << endl;
		cout << " " << endl;
		cout << "OffDiag_C0C1: " << offdiagC0_C1(0,0) << " " << offdiagC0_C1(0,1) << " " << offdiagC0_C1(0,2) << endl;
		cout << "OffDiag_C0C1: " << offdiagC0_C1(1,0) << " " << offdiagC0_C1(1,1) << " " << offdiagC0_C1(1,2) << endl;
		cout << "OffDiag_C0C1: " << offdiagC0_C1(2,0) << " " << offdiagC0_C1(2,1) << " " << offdiagC0_C1(2,2) << endl;
		cout << " " << endl;
		cout << "Diag_C0: " << (diag[c0])(0,0) << " " << (diag[c0])(0,1) << " " << (diag[c0])(0,2) << endl;
		cout << "Diag_C0: " << (diag[c0])(1,0) << " " << (diag[c0])(1,1) << " " << (diag[c0])(1,2) << endl;
		cout << "Diag_C0: " << (diag[c0])(2,0) << " " << (diag[c0])(2,1) << " " << (diag[c0])(2,2) << endl;
		cout << " " << endl;
		cout << "Diag_C1: " << (diag[c1])(0,0) << " " << (diag[c1])(0,1) << " " << (diag[c1])(0,2) << endl;
		cout << "Diag_C1: " << (diag[c1])(1,0) << " " << (diag[c1])(1,1) << " " << (diag[c1])(1,2) << endl;
		cout << "Diag_C1: " << (diag[c1])(2,0) << " " << (diag[c1])(2,1) << " " << (diag[c1])(2,2) << endl;
		cout << " " << endl;
		cout << "Residual_C0: " << (rCell[c0])[0] << endl;
		cout << "Residual_C0: " << (rCell[c0])[1] << endl;
		cout << "Residual_C0: " << (rCell[c0])[2] << endl;
		cout << " " << endl;
		cout << "Residual_C1: " << (rCell[c1])[0] << endl;
		cout << "Residual_C1: " << (rCell[c1])[1] << endl;
		cout << "Residual_C1: " << (rCell[c1])[2] << endl;
		cout << " " << endl;*/
		
	      }
	  }

      }

  }
 private:
  
  const GeomFields& _geomFields;
  Field& _varField;
  const T_Scalar _RRConstant;
  const T_Scalar _interfaceUnderRelax;
  const bool _Anode;
  const bool _Cathode;
  const bool _bInterfaceHeatSource;
  
};

#endif
