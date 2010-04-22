#ifndef _TUNNELINGDISCRETIZATION_H_
#define _TUNNELINGDISCRETIZATION_H_

#include "PhysicsConstant.h"
#include "ElectricBC.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "CRMatrix.h"
#include "DielectricOneDimColumn.h"
#include "ElectricUtilityFunctions.h"

//Diag type: 2x2Tensor
/***  | d0,  d1 |
      | d2,  d3 |    ***/
//OffDiag type: 2x2Tensor
/***  | o0,  o1 |
      | o2,  o3 |    ***/
//X type: VectorT2
/***  | x0 |
      | x1 |         ***/

// in tunneling model, only d0 and x0 are modified
// the rest remains unchanged

template <class X, class Diag, class OffDiag>
class TunnelingDiscretization : public Discretization
{

 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;

  TunnelingDiscretization(const MeshList& meshes,
			  const GeomFields& geomFields,
			  const Field& varField,
			  const Field& varN1Field,
			  const Field& totaltrapsField,
			  const Field& conductionbandField,
			  const ElectricModelConstants<T_Scalar>& constants,
			  Field& transmission, 
			  vector<shared_ptr<DielectricOneDimColumn<T_Scalar> > >& columnList) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varN1Field(varN1Field),
    _totaltrapsField(totaltrapsField),
    _conductionbandField(conductionbandField),
    _constants(constants),
    _transmission(transmission),
    _columnList(columnList)
      {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    
    const StorageSite& cells = mesh.getCells();
    
    const int nCells = cells.getSelfCount();
    
    const TArray& electron_totaltraps = dynamic_cast<const TArray&> (_totaltrapsField[cells]);

    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    TArray& transmission = dynamic_cast<TArray&> (_transmission[cells]);
   
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(_varField[cells]);
    
    const XArray& xN1Cell = dynamic_cast<const XArray&>(_varN1Field[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    OffDiagArray& offdiag = matrix.getOffDiag();
    
    T_Scalar energystep = 0.01;
    T_Scalar currentTime = 1.0;
    double fluxCoeff, fermilevel, scatterfactor;

    //const T_Scalar volt = getMembraneVoltage(currentTime);
    const T_Scalar volt = 150;

    const T_Scalar& electron_effmass = _constants["electron_effmass"];
    const T_Scalar& temperature = _constants["OP_temperature"];
    const T_Scalar& electron_capture_cross = _constants["electron_capture_cross"];
    const T_Scalar& electron_trapdepth = _constants["electron_trapdepth"];
    const T_Scalar& fermilevelsubstrate = -_constants["substrate_workfunction"];
    const T_Scalar& fermilevelmembrane = -_constants["substrate_workfunction"] - volt;
    const T_Scalar& dielectric_ionization = _constants["dielectric_ionization"];

    const T_Scalar alpha = 4.0 * PI * (electron_effmass*ME) / pow(H_SI, 3.0);

    //for the tunneling from substrate to traps and from traps to substrate, the same fermilevel is used
    //thus the same transmission can be used twice
    //the flux coeff is a little different

    fermilevel = fermilevelsubstrate;

    
    for (T_Scalar en = fermilevel-1.0; en <= fermilevel+1.0; en += energystep){

      const T_Scalar supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);

      //cout << "supplyfunction" << supplyfunction << endl;

      const T_Scalar fermifunction = FermiFunction(en, fermilevel, temperature);

      //cout << "fermifunction" << fermifunction << endl;

      const string flag = "substrate";

      ElectronTransmissionCoefficient(en, transmission, conduction_band,
      			      dielectric_ionization, electron_effmass, flag, _columnList);
      

      for(int c=0; c<nCells; c++){
	/*
	//debug transmission interpolation
	T p = (cellCentriod[c][2] - 0.0) / 0.1e06 * 150.0;
	T deltaX = 0.1e-6 / 20.0;
	T cond = -(p + 3.0);
	T factor = -2.0/HBAR * sqrt(2.0*electron_effmass*ME*QE);
	T exponent = factor * sqrt(PositiveValueOf(cond - en)) * deltaX;
	*/
	transmission[c] = 1.0;
	
	const T_Scalar stcap = electron_capture_cross * cellVolume[c]; 

	const T_Scalar endiff = en - (conduction_band[c]-electron_trapdepth);
	
	//tunneling from substrate to traps
	if (endiff < 0)
	  scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	else scatterfactor = 1.0;

	//cout << "scatter " << scatterfactor << endl;

	fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * fermifunction * scatterfactor * energystep * QE;

	//r = coeff * (N-n) explicit; no diag;
	//if implicit, the xN1 is needed, and a fraction factor

	rCell[c][0] += (fluxCoeff * (electron_totaltraps[c] - xN1Cell[c][0])); 
		
	//diag[c][0] += fluxCoeff;

	//tunneling from traps to substrate
	if (endiff > 0)
	  scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	else scatterfactor = 1.0;

	fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * (1-fermifunction) * scatterfactor * energystep * QE;
	
	//r = coeff * (-n) explicit; no diag;
	//rCell[c][0] += (fluxCoeff * (- xN1Cell[c][0])); 

	//diag[c][0] += fluxCoeff;
      }
    }
    
    
     //for the tunneling from membrane to traps and from traps to membrane, the same fermilevel is used
    //thus the same transmission can be used twice
    //the flux coeff is a little different
    /*
    fermilevel = fermilevelmembrane;
    
    for (T_Scalar en = fermilevel-4.0; en <= fermilevel+4.0; en += energystep){

      const T_Scalar supplyfunction = ElectronSupplyFunction(en, fermilevel, temperature);

      const T_Scalar fermifunction = FermiFunction(en, fermilevel, temperature);

      const string flag = "membrane";

      ElectronTransmissionCoefficient(en, transmission, conduction_band,
				      dielectric_ionization, electron_effmass, flag, _columnList);
      
      for(int c=0; c<nCells; c++){
	//transmission[c]=1.0;
	
	const T_Scalar stcap = electron_capture_cross * cellVolume[c];  

	const T_Scalar endiff = en - (conduction_band[c]-electron_trapdepth);
	
	//tunneling from membrane to traps
	
	if (endiff < 0)
	  scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	else scatterfactor = 1.0;

	fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * fermifunction * scatterfactor * energystep * QE;
	if (c==0) cout << fluxCoeff << endl;	
	rCell[c][0] += (fluxCoeff * (electron_totaltraps[c] - xN1Cell[c][0])); 

	//diag[c] += fluxCoeff;

	//tunneling from traps to membrane

	if (endiff > 0)
	  scatterfactor = exp(-QE * fabs(endiff)/(K_SI*temperature));
	else scatterfactor = 1.0;
	
	fluxCoeff = alpha * stcap * transmission[c] * supplyfunction * (1-fermifunction) * scatterfactor * energystep * QE;

	rCell[c][0] += (fluxCoeff * (-xN1Cell[c][0])); 

	//diag[c][0] += fluxCoeff;
	}
    }
  
    */
  }
 
 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _totaltrapsField;
  const Field& _conductionbandField;
  const ElectricModelConstants<T_Scalar>& _constants;
  Field& _transmission;
  vector<shared_ptr<DielectricOneDimColumn<T_Scalar> > >& _columnList;
};


#endif
