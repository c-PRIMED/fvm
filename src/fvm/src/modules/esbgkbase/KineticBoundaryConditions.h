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
  
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;

  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef FluxJacobianMatrix<Diag,X> FMatrix;
  typedef DiagonalMatrix<Diag,X> BBMatrix;

  typedef Array<Diag> DiagArray;
  typedef Array<OffDiag> OffDiagArray;

  typedef Array<X> XArray;
  typedef Array<VectorT3> VectorT3Array;


  KineticBoundaryConditions(const StorageSite& faces,
			    const Mesh& mesh,
			    const GeomFields& geomFields,
			    const Quadrature<X>& quadrature,
			    DistFunctFields<X>& dsfPtr,
			    Field& varField,
			    Field& fluxField,
			    MultiFieldMatrix& matrix,
			    MultiField& xField,
			    MultiField& rField,
			    int direction) :
  _faces(faces),
    _cells(mesh.getCells()),
    _ibType(mesh.getIBType()),
    _quadrature(quadrature),
    _dsfPtr(dsfPtr),
    _faceCells(mesh.getFaceCells(_faces)),
    _varField(varField),
    _fluxField(fluxField),
    _xIndex(&_varField,&_cells),
    _fluxIndex(&_fluxField,&_faces),
    _dRdX(dynamic_cast<CCMatrix&>(matrix.getMatrix(_xIndex,_xIndex))),
    _dFluxdX(dynamic_cast<FMatrix&>(matrix.getMatrix(_fluxIndex,_xIndex))),
    _dFluxdFlux(dynamic_cast<BBMatrix&>(matrix.getMatrix(_fluxIndex,_fluxIndex))),
    _assembler(_dRdX.getPairWiseAssembler(_faceCells)),
    _dRdXDiag(_dRdX.getDiag()),
    _x(dynamic_cast<XArray&>(xField[_xIndex])),
    _r(dynamic_cast<XArray&>(rField[_xIndex])),
    _flux(dynamic_cast<XArray&>(xField[_fluxIndex])),
    _rFlux(dynamic_cast<XArray&>(rField[_fluxIndex])),
    _areaMagField(geomFields.areaMag),
    _faceAreaMag(dynamic_cast<const TArray&>(_areaMagField[_faces])),
    _areaField(geomFields.area),
    _faceArea(dynamic_cast<const VectorT3Array&>(_areaField[_faces])),
    _direction(direction)
      {}

  void applyWallBC(int f,const X& WallVelocity,const X& WallTemperature) const
  {

    const double pi(acos(0.0)); 
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,0);

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
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	const X fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	if (c_dot_en > T_Scalar(0.0))
	  {
	    Nmr = Nmr + fnd[c0];
	  }
	else
	  {
	    Dmr = Dmr + fwall;
	  }	
      }
    const X nwall = Nmr/Dmr; // wall number density for initializing Maxwellian

    const X fwall = nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[_direction]-uwall,2.0)+pow(cy[_direction]-vwall,2.0)+pow(cz[_direction]-wwall,2.0))/Twall);

    // the current value of flux and its Jacobians
    const X fluxB = -_r[c1];
    const OffDiag dFluxdXC0 = -_assembler.getCoeff10(f);
    const Diag dFluxdXC1 = -_dRdXDiag[c1];
    const OffDiag dRC0dXC1 = _assembler.getCoeff01(f);

    // since we know the boundary value, compute the boundary
    // x correction and it's contribution to the residual for c0;
    // we can then eliminate the coefficient to the boundary cell

    const X dXC1 = fwall - _x[c1];
    const X dFlux = dFluxdXC1*dXC1;
    const X dRC0 = dRC0dXC1*dXC1;
    _r[c0] += dRC0;

    _assembler.getCoeff01(f) = NumTypeTraits<OffDiag>::getZero();

    // set the boundary value and make its correction equation an
    // identity
    _x[c1] = fwall;
    _assembler.getCoeff10(f) = NumTypeTraits<OffDiag>::getZero();
    _r[c1] = NumTypeTraits<X>::getZero();
    _dRdXDiag[c1] = NumTypeTraits<Diag>::getNegativeUnity();

    // setup the equation for the boundary flux correction
    _dFluxdX.setCoeffL(f,dFluxdXC0);
    _dFluxdX.setCoeffR(f,NumTypeTraits<OffDiag>::getZero());
    _flux[f] = fluxB;
    _rFlux[f] = dFlux;
    _dFluxdFlux[f] = NumTypeTraits<Diag>::getNegativeUnity();
  
  }

  //  void applySpecularBC(int f) const
  //{
  //}

 protected:
  const StorageSite& _faces;
  const StorageSite& _cells;
  const Array<int>& _ibType;
  const Quadrature<X>& _quadrature;
  const DistFunctFields<X>& _dsfPtr;
  const CRConnectivity& _faceCells;
  const Field& _varField;
  const Field& _fluxField;
  const MultiField::ArrayIndex _xIndex;
  const MultiField::ArrayIndex _fluxIndex;
  CCMatrix& _dRdX;
  FMatrix& _dFluxdX;
  BBMatrix& _dFluxdFlux;
  CCAssembler& _assembler;
  DiagArray& _dRdXDiag;
  XArray& _x;
  XArray& _r;
  XArray& _flux;
  XArray& _rFlux;
  const Field& _areaMagField;
  const TArray& _faceAreaMag;
  const Field& _areaField;
  const VectorT3Array& _faceArea;
  const int _direction;
};
  
#endif
