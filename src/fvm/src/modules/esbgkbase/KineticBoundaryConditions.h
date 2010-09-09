#ifndef _KINETICBOUNDARYCONDITIONS_H_
#define _KINETICBOUNDARYCONDITIONS_H_

#include "Mesh.h"

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

template<class X,class Diag,class OffDiag>
class KineticBoundaryConditions
{
 public :
  
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;

  typedef Array<int> IntArray;
  
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;

  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef FluxJacobianMatrix<Diag,X> FMatrix;
  typedef DiagonalMatrix<Diag,X> BBMatrix;

  typedef Array<Diag> DiagArray;
  typedef Array<OffDiag> OffDiagArray;

  typedef Array<X> XArray;
  typedef Vector<X,3> VectorX3; //new for esbgk

  typedef Array<VectorT3> VectorT3Array;


  KineticBoundaryConditions(const StorageSite& faces,
			    const Mesh& mesh,
			    const GeomFields& geomFields,
			    const Quadrature<X>& quadrature,
			    DistFunctFields<X>& dsfPtr):
			  
    _faces(faces),
    _cells(mesh.getCells()),
    _ibType(dynamic_cast<const IntArray&>(geomFields.ibType[_cells])),
    _quadrature(quadrature),
    _dsfPtr(dsfPtr),
    _faceCells(mesh.getFaceCells(_faces)),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces]))
    
  {}
    
  void applyDiffuseWallBC(int i,const VectorX3&  WallVelocity,const X& WallTemperature) const
  {
    
    const double pi(acos(0.0)); 
    const int c0 = _faceCells(i,0);
    const int c1 = _faceCells(i,0);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    // For the cell c0, loop over all directions
    // check for c.n sign and add to numerator part
    // or denominator part (to compute wall number
    // density). The wall distribution function (through the 
    // wall number density and wall temperature) acts like
    // the variable bValue in the fvmbase module
    const int numDirections = _quadrature.getDirCount();
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    
    const X uwall = WallVelocity[0];
    const X vwall = WallVelocity[1];
    const X wwall = WallVelocity[2];
    
    const X Twall = WallTemperature;
    
    X Nmr(0.0) ;
    X Dmr(0.0) ;
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& f = dynamic_cast< TArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[i]/_faceAreaMag[i];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	const X fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	if (c_dot_en > T_Scalar(0.0))
	  {
	    Nmr = Nmr + f[c0];
	  }
	else
	  {
	    Dmr = Dmr + fwall;
	  }	
      }
    const X nwall = Nmr/Dmr; // wall number density for initializing Maxwellian
    
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& f = dynamic_cast< TArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[i]/_faceAreaMag[i];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	if (c_dot_en < T_Scalar(0.0))
	  {
	    f[c0] = nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	    
	  }
	
      }
  }

void applyDiffuseWallBC(const VectorX3& bVelocity,const X& bTemperature) const
  {
    for (int i=0; i<_faces.getCount();i++)
    applyDiffuseWallBC(i,bVelocity,bTemperature);
  }
  
  
  void applyDiffuseWallBC(const FloatValEvaluator<VectorX3>& bVelocity,const FloatValEvaluator<X>& bTemperature) const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyDiffuseWallBC(i,bVelocity[i],bTemperature[i]);

  }

 
  
 protected:

  const StorageSite& _faces;
  const StorageSite& _cells;
  const Array<int>& _ibType;
  const Quadrature<X>& _quadrature;
  const DistFunctFields<X>& _dsfPtr;
  const CRConnectivity& _faceCells;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  //const int _direction;
};
  
#endif
