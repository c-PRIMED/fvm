#ifndef _INJECTIONDISCRETIZATION_H_
#define _INJECTIONDISCRETIZATION_H_

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

/**************************************************
Diag type: 2x2Tensor
      | d00,  d01 |
      | d10,  d11 |    

OffDiag type: 2x2Tensor
      | o00,  o01 |
      | o10,  o11 |    

X type: VectorT2
      | x0 |
      | x1 |         

"0" is trapped charge
"1" is band charge

Injection model only modifies x1 

*************************************************/

template <class X, class Diag, class OffDiag>
class InjectionDiscretization : public Discretization
{

 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;
  typedef Array<Vector<T_Scalar, 3> > VectorT3Array;


  InjectionDiscretization(const MeshList& meshes,
			  const GeomFields& geomFields,
			  const Field& varField,
			  const Field& electricField,
			  const Field& conductionbandField,
			  const ElectricModelConstants<T_Scalar>& constants
			  ) :
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _electricField(electricField),
    _conductionbandField(conductionbandField),
    _constants(constants)    
      {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const int nCells = cells.getSelfCount();

    const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array&> (_electricField[cells]);
    
    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    //TArray& transmission = dynamic_cast<TArray&> (_transmissionField[cells]);

    shared_ptr<TArray> ts(new TArray(nCells));
    *ts = 0;
    TArray& transmission = *ts;

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
     
    const VectorT3Array& cellCentroid = 
      dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);
    
    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    OffDiagArray& offdiag = matrix.getOffDiag();

    
    const T_Scalar& electron_effmass = _constants["electron_effmass"];
    const T_Scalar& temperature = _constants["OP_temperature"];
    const T_Scalar& fermilevelsubstrate = -_constants["substrate_workfunction"];
    //const T_Scalar& fermilevelmembrane = -_constants["substrate_workfunction"] - volt;
    
    const int& subID = _constants["substrate_id"];
    const int& memID = _constants["membrane_id"];
    const int& nLevel = _constants["nLevel"];
    const int& normal = _constants["normal_direction"];

    T_Scalar fluxCoeff(0), fermilevel(0), scatterfactor(0);
    
    const T_Scalar epsilon = 1e-18;
    
    const T_Scalar alpha = 4.0 * PI * (electron_effmass*ME) / pow(H_SI, 3.0);
    
    for(int c=0; c<nCells; c++){
      transmission[c] = 0.0;
    } 

    fermilevel = fermilevelsubstrate;   

    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;
	const StorageSite& faces = fg.site;
	
	if (fg.id == subID){
	  
	  const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	  const CRConnectivity& cellCells = mesh.getCellCells();
	  const int nFaces = faces.getCount();
	  
	  for(int f=0; f<nFaces; f++){

	    int c0 = faceCells(f,0);
	    int c1 = faceCells(f,1);
	    	     
	    (transmission)[c1] = 1.0;

	    int low = c1;
	    int me = c0;	      
	    int high = c0;	       
	    
	    //for(int l=0; l<nLevel; l++){
	    for(int l=0; l< int(nLevel/10); l++){

	      T_Scalar enMe = conduction_band[me];
	      T_Scalar enLow = conduction_band[low];
	      T_Scalar deltaEn = fabs(enMe - enLow);
	      T_Scalar dX = fabs(cellCentroid[me][normal] - cellCentroid[low][normal]);

	      // cout << fermilevel << "  " << enMe << endl;
	     
	      //calcualte transmission[me] from start point to me
	      //--------------------------------------------------------------------------------
	      int in_low = c1;
	      int in_me = c0;	      
	      int in_high = c0;
	      int count = 0;

	      T_Scalar deltaX = cellCentroid[in_me][normal] - cellCentroid[in_low][normal];
	      T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	      T_Scalar valueMe = PositiveValueOf( conduction_band[in_me] - enMe);
	      T_Scalar valueLow = PositiveValueOf( conduction_band[in_low] - enMe);
	      T_Scalar avg = (valueMe + valueLow) / 2.0;
	      T_Scalar exponent = factor * sqrt(avg) * fabs(deltaX);
		
	      (transmission)[in_me] = (transmission)[in_low] * exp(exponent);
	      

	      while (in_me != me && count <= l){
		
		const int nbc = cellCells.getCount(in_me);

		for(int nc = 0; nc < nbc; nc++){

		  const int neighbor = cellCells(in_me, nc);
		  
		  const T_Scalar dr = cellCentroid[in_me][normal] - cellCentroid[neighbor][normal];
		 
		  if ((fabs(dr) > epsilon) && (dr < 0.0) && (neighbor < nCells)) {

		    in_high = neighbor;
		    in_low = in_me;
		    in_me = in_high;
		  }
		}
	
		T_Scalar deltaX = cellCentroid[in_me][normal] - cellCentroid[in_low][normal];
		T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
		T_Scalar valueMe = PositiveValueOf( conduction_band[in_me] - enMe);
		T_Scalar valueLow = PositiveValueOf( conduction_band[in_low] - enMe);
		T_Scalar avg = (valueMe + valueLow) / 2.0;
		T_Scalar exponent = factor * sqrt(avg) * fabs(deltaX);
		
		(transmission)[in_me] = (transmission)[in_low] * exp(exponent);
		
		count ++;
	      }
	      //--------------------------------------------------------------------------------


	      if (in_me != me)
		throw CException ("wrong injection cell loop!");

	      const T_Scalar supplyfunction = ElectronSupplyFunction(enMe, fermilevel, temperature);

	      const T_Scalar fermifunction = FermiFunction(enMe, fermilevel, temperature);
		
     	      fluxCoeff = alpha * (transmission)[me] * supplyfunction * fermifunction * fabs(deltaEn) * QE * cellVolume[me] / fabs(dX) ;
	      
	      //fluxCoeff = alpha * (transmission)[me] * supplyfunction * fermifunction * fabs(electric_field[me][normal]) * cellVolume[me] * QE;

	      fluxCoeff = fabs(fluxCoeff);
	      
	      //cout << me <<" "<< (*transmission)[me] << "  " << supplyfunction << "  " << fermifunction << "   " << fluxCoeff << endl;	      

	      rCell[me][1] += fluxCoeff; 
	      
	      const int nbc = cellCells.getCount(me);
	      
	      for(int nc = 0; nc < nbc; nc++){

		const int neighbor = cellCells(me, nc);
		  
		const T_Scalar dr = cellCentroid[me][normal] - cellCentroid[neighbor][normal];
		 
		if ((fabs(dr) > epsilon) && (dr < 0.0) && (neighbor < nCells)) {
		  high = neighbor;
		  low = me;
		  me = high;
		}
	      }
	    }
	  }	
	}
      }
  }
	


 private:
  const GeomFields& _geomFields;
  const Field& _conductionbandField;
  const Field& _varField;
  const Field& _electricField;
  const ElectricModelConstants<T_Scalar>& _constants;
  

};





#endif
