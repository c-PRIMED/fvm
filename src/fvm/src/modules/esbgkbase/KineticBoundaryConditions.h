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
    
  void applyDiffuseWallBC(int f,const VectorX3&  WallVelocity,const X& WallTemperature) const
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
    //cout << "uwall " << uwall << endl;
    const X Twall = WallTemperature;
    
    X Nmr(0.0) ;
    X Dmr(0.0) ;
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	const X fwall = 1.0/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	if (c_dot_en -wallV_dot_en > T_Scalar(0.0))
	  {
	    Nmr = Nmr + dsf[c0];
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
	TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];	
	const X wallV_dot_en = uwall*en[0]+vwall*en[1]+wwall*en[2];
	if (c_dot_en-wallV_dot_en < T_Scalar(0.0))
	  {
	    dsf[c1] = nwall/pow(pi*Twall,1.5)*exp(-(pow(cx[j]-uwall,2.0)+pow(cy[j]-vwall,2.0)+pow(cz[j]-wwall,2.0))/Twall);
	    
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

  void applySpecularWallBC(int f) const
  {
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,0);

    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    // For the cell c0, loop over all directions
    // check for c.n sign and compute the incident molecules
    // velocity direction if the current direction is the reflected 
    // molecules direction
    const int numDirections = _quadrature.getDirCount();
    const TArray& cx = dynamic_cast<const TArray&>(*_quadrature.cxPtr);
    const TArray& cy = dynamic_cast<const TArray&>(*_quadrature.cyPtr);
    const TArray& cz = dynamic_cast<const TArray&>(*_quadrature.czPtr);
    const int N1 = _quadrature.getNVCount();
    const int N2 = _quadrature.getNthetaCount();
    const int N3 = _quadrature.getNphiCount();
    const X dcx = (cx[N1-1]-cx[0])/(N1-1);
    const X dcy = (cy[N1-1]-cy[0])/(N2-1);
    const X dcz = (cz[N1-1]-cz[0])/(N3-1);
    for (int j=0; j<numDirections; j++)
      {
	// Find incident molecule directions
	const VectorT3 en = _faceArea[f]/_faceAreaMag[f];
	const X  c_dot_en = cx[j]*en[0]+cy[j]*en[1]+cz[j]*en[2];
	
	
	if(c_dot_en < T_Scalar(0.0)) //incoming molecules to interior
	  {
	    const X cx_incident = cx[j] - 2.0*c_dot_en*en[0];
	    const X cy_incident = cy[j] - 2.0*c_dot_en*en[1];
	    const X cz_incident = cz[j] - 2.0*c_dot_en*en[2];
	    
	    // Find incident molecule indices
	    // Right now only for the case of uniform cartesian velocity grid
	    // A search has to be performed to generalize it to the non-uniform
	    // velocity grids where dcx, dcy, dcz are not constant
	   
	    	   
	    
	    const int i_incident = (cx_incident-cx[0])/dcx;
	    const int j_incident = (cy_incident-cy[0])/dcy;
	    const int k_incident = (cz_incident-cz[0])/dcz;
	    const int direction_incident = k_incident+N3*j_incident+N3*N2+i_incident;
	
	    // The value of 'f' in the incident direction will serve
	    // as the boundary value for the current(reflected) direction
	    Field& fndw = *_dsfPtr.dsf[j];
	    TArray& dsfw = dynamic_cast<TArray&>(fndw[_cells]);
	    
	    Field& fnd = *_dsfPtr.dsf[direction_incident];
	    const TArray& dsf = dynamic_cast<const TArray&>(fnd[_cells]);
	    
	    dsfw[c1] = dsf[c0]; //write into boundary
	  }
      }
  }   
  void applySpecularWallBC() const
  {
    for (int i=0; i<_faces.getCount();i++)
      applySpecularWallBC(i);
  }

 
  void applyCopyWallBC(int f) const
  {
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,0);
    
    if (_ibType[c0] != Mesh::IBTYPE_FLUID)
      return;
    
    const int numDirections = _quadrature.getDirCount();
    
    for (int j=0; j<numDirections; j++)
      {
	Field& fnd = *_dsfPtr.dsf[j];
	TArray& dsf = dynamic_cast< TArray&>(fnd[_cells]);
	dsf[c1]=dsf[c0];
      }
  }
  void applyCopyWallBC() const
  {
    for (int i=0; i<_faces.getCount();i++)
      applyCopyWallBC(i);
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
  
};
  
#endif
