#ifndef _TRAPBANDTUNNELINGDISCRETIZATION_H_
#define _TRAPBANDTUNNELINGDISCRETIZATION_H_

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

template <class X, class Diag, class OffDiag>
class TrapBandTunnelingDiscretization : public Discretization
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef Array<T_Scalar> TArray;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;
  typedef Array<X> XArray;
  typedef Array<Vector<T_Scalar, 3> > VectorT3Array;

  TrapBandTunnelingDiscretization(const MeshList& meshes,
				  const GeomFields& geomFields,
				  const Field& varField,
				  const Field& varN1Field,
				  const Field& electricField,
				  const Field& conductionbandField,
				  const ElectricModelConstants<T_Scalar>& constants):
    Discretization(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _varN1Field(varN1Field),
    _electricField(electricField),
    _conductionbandField(conductionbandField),
    _constants(constants)
      {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const TArray& conduction_band = dynamic_cast<const TArray&> (_conductionbandField[cells]);

    const VectorT3Array& electric_field = dynamic_cast<const VectorT3Array&> (_electricField[cells]);

    const VectorT3Array& cellCentroid = 
      dynamic_cast<const VectorT3Array& > (_geomFields.coordinate[cells]);

    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,cVarIndex));
    
    const XArray& xCell = dynamic_cast<const XArray&>(_varField[cells]);
    
    const XArray& xN1Cell = dynamic_cast<const XArray&>(_varN1Field[cells]);
    
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);
       
    DiagArray& diag = matrix.getDiag();

    OffDiagArray& offdiag = matrix.getOffDiag();

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);

    const T_Scalar& electron_trapdepth = _constants["electron_trapdepth"];
    const T_Scalar& electron_capture_cross = _constants["electron_capture_cross"];
    const T_Scalar& electron_effmass = _constants["electron_effmass"];

    const int subID = 5;
    const int memID = 4;
    const int nLevel = 100;


    //start from substrate boundary faces
    foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
    {
      const FaceGroup& fg = *fgPtr;

      if (fg.id == subID) {

	const StorageSite& faces = fg.site;
	const CRConnectivity& faceCells = mesh.getFaceCells(faces);

	for(int f=0; f<faces.getCount(); f++){
	  //for (int f=0; f<1; f++){
	  int c0 = faceCells(f,0);
	  int c1 = faceCells(f,1);
	   
	  int low = c1;
	  int me = c0;	      
	  int high = c0 - 4;	

	  for(int l=0; l<nLevel; l++){	
	    T_Scalar en = conduction_band[me] - electron_trapdepth;
	    T_Scalar transmissionLow = 1;
	    T_Scalar transmissionHigh = 1;
	    int jLow = l,  jHigh = l;            //level index
	    int iLow = low, iHigh = high;       // real cell index
	    int lowFound = 0, highFound = 0;    // found or not
	   
	    T_Scalar ef = mag(electric_field[me]);
	    T_Scalar factor = cellVolume[me] * QE * ef * ef * electron_capture_cross / 
	      (16 * PI*PI * HBAR_SI * electron_effmass * electron_trapdepth);
	    //checking conduction band edge on membrane
	    jHigh = l+1;
	    while(jHigh < nLevel){
	      T_Scalar delta = cellCentroid[iHigh][2] - cellCentroid[me][2];
	      T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	      T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
	      T_Scalar valueHigh = PositiveValueOf( conduction_band[iHigh] - en);
	      T_Scalar avg = (valueMe + valueHigh) / 2.0;
	      T_Scalar exponent = factor * sqrt(avg) * delta; 
	      transmissionHigh *= exp(exponent);
	      if(en-conduction_band[iHigh] >0){	
		highFound = 1;
		//cout << jHigh << "  " << iHigh << "  " <<transmissionHigh << endl;
		break;
	      }
	      jHigh ++;
	      iHigh = iHigh - 4;
	    }
	    jLow = l-1;
	    while (jLow >= 0 ){
	      T_Scalar delta = cellCentroid[me][2] - cellCentroid[iLow][2];
	      T_Scalar factor = -2.0/HBAR_SI * sqrt(2.0*electron_effmass*ME*QE);
	      T_Scalar valueMe = PositiveValueOf( conduction_band[me] - en);
	      T_Scalar valueLow = PositiveValueOf( conduction_band[iLow] - en);
	      T_Scalar avg = (valueMe + valueLow) / 2.0;
	      T_Scalar exponent = factor * sqrt(avg) * delta; 
	      transmissionLow *= exp(exponent);
	      if(en-conduction_band[iLow] >0){	
		lowFound = 1;
		//cout << jLow << "  " << iLow << "  " <<transmissionLow << endl;
		break;
	      }
	      jLow --;
	      iLow = iLow + 4;
	    }
	    if( lowFound == 1 || highFound == 1){
	      rCell[me][0] -= factor * (transmissionLow + transmissionHigh) * xCell[me][0];
	      diag[me][0] -= factor  * (transmissionLow + transmissionHigh);
	    }
	    
	    if (lowFound == 1){
	      rCell[iLow][1] += factor * transmissionLow * xCell[me][0];
	      diag[iLow][1] += factor * transmissionLow;
	    }

	    if (highFound == 1){
	      rCell[iHigh][1] += factor * transmissionHigh * xCell[me][0];
	      diag[iHigh][1] += factor * transmissionHigh;
	    }
	    
	    low = me;
	    me = high;
	    high = high - 4;
	    
	  }
	}
      }
    }
  }

 private:
  const GeomFields& _geomFields;
  const Field& _varField;
  const Field& _varN1Field;
  const Field& _electricField;
  const Field& _conductionbandField;
  const ElectricModelConstants<T_Scalar>& _constants;
 
	      
};

#endif
