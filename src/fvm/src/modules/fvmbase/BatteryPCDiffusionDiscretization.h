// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BATTERYPCDIFFUSIONDISCRETIZATION_H_
#define _BATTERYPCDIFFUSIONDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "Gradient.h"
#include "DiagonalTensor.h"

template<class T>
inline T harmonicAverageVector(const T& x0, const T& x1)
{
  const T sum = x0+x1;
  T harmAvg = sum;
  int vectorLength = 0;
  NumTypeTraits<T>::getShape(&vectorLength);
  for (int v=0; v<vectorLength; v++)
    {
      if (sum[v] != 0.0)
	harmAvg[v] = 2.0*x0[v]*x1[v]/sum[v];
    }
  return harmAvg;
}

template<class X, class Diag, class OffDiag>
class BatteryPCDiffusionDiscretization : public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::PairWiseAssembler CCAssembler;

  typedef Gradient<X> XGrad;
  //typedef Array<int> IntArray;
  
  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<XGrad> GradArray;

  ///// ASSUMING ONLY ONE SPECIES, BELOW IS TRUE //////
  // X = VectorT2
  // XArray = VectorT2Array;
  // XGrad = Gradient<VectorT2>
  // GradArray = Array<Gradient<VectorT2>>

  ///// ASSUMING ONLY ONE SPECIES AND TEMP INCLUDED, BELOW IS TRUE //////
  // X = VectorT3
  // XArray = VectorT3Array;
  // XGrad = Gradient<VectorT3>
  // GradArray = Array<Gradient<VectorT3>>

  typedef Gradient<T_Scalar> TGrad;
  
  BatteryPCDiffusionDiscretization(const MeshList& meshes,
                          const GeomFields& geomFields,
                          Field& varField,
                          const Field& diffusivityField,
                          const Field& varGradientField) :
      Discretization(meshes),
      _geomFields(geomFields),
      _varField(varField),
      _diffusivityField(diffusivityField),
      _varGradientField(varGradientField)
  {}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();
    
    const MultiField::ArrayIndex cVarIndex(&_varField,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(cVarIndex,
                                                             cVarIndex));    

    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[cells]);
   
    const TArray& cellVolume =
      dynamic_cast<const TArray&>(_geomFields.volume[cells]);
  
    DiagArray& diag = matrix.getDiag();

    const XArray& xCell = dynamic_cast<const XArray&>(xField[cVarIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[cVarIndex]);

    const GradArray& xGradCell =
      dynamic_cast<const GradArray&>(_varGradientField[cells]);

    const XArray& diffCell =
      dynamic_cast<const XArray&>(_diffusivityField[cells]);
   
    foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
      {
	const FaceGroup& fg = *fgPtr;
	const StorageSite& faces = fg.site;
	const int nFaces = faces.getCount();
	const CRConnectivity& faceCells = mesh.getFaceCells(faces);
	const VectorT3Array& faceArea =
	  dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);    
	const TArray& faceAreaMag =
	  dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
	//const VectorT3Array& faceCentroid = dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[faces]);
	CCAssembler& assembler = matrix.getPairWiseAssembler(faceCells);
	
	int XLength = 0; 
	NumTypeTraits<X>::getShape(&XLength);

	for(int f=0; f<nFaces; f++)
	  {
	    const int c0 = faceCells(f,0);
	    const int c1 = faceCells(f,1);

	    T_Scalar vol0 = cellVolume[c0];
	    T_Scalar vol1 = cellVolume[c1];
        
	    VectorT3 ds=cellCentroid[c1]-cellCentroid[c0];
    

	    // for ib faces ignore the solid cell and use the face centroid for diff metric
	    /*
	      if (((ibType[c0] == Mesh::IBTYPE_FLUID)
	      && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
	      ((ibType[c1] == Mesh::IBTYPE_FLUID)
	      && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
	      {
	      if (ibType[c0] == Mesh::IBTYPE_FLUID)
	      {
	      vol1 = 0.;
	      ds = faceCentroid[f]-cellCentroid[c0];
	      }
	      else
	      {
	      vol0 = 0.;
	      ds = cellCentroid[c1]-faceCentroid[f];
	      }
	      }*/

	    //for(int v=0; v<XLength; v++)
	      
	    X faceDiffusivity(NumTypeTraits<X>::getZero());
	    if (vol0 == 0.)
	      faceDiffusivity = (diffCell[c1]);
	    else if (vol1 == 0.)
	      faceDiffusivity = (diffCell[c0]);
	    else
	      {
		faceDiffusivity = harmonicAverageVector(diffCell[c0],diffCell[c1]);
	      }
		
	    const T_Scalar diffMetric = faceAreaMag[f]*faceAreaMag[f]/dot(faceArea[f],ds);

	    //do things element-wise manually
	    // Each 'v' below is on equation
	    // first equation is potential, followed by species equation, then thermal
	    for(int v=0; v<XLength; v++)
	      {
		const T_Scalar diffCoeff = (faceDiffusivity[v])*diffMetric;

		const VectorT3 secondaryCoeff = (faceDiffusivity[v])*(faceArea[f]-ds*diffMetric);
        
		//extract this v's gradient
		const XGrad gradF_X = ((xGradCell[c0]*vol0+xGradCell[c1]*vol1)/(vol0+vol1));
		TGrad gradF(NumTypeTraits<TGrad>::getZero());
		gradF[0] = gradF_X[0][v];
		gradF[1] = gradF_X[1][v];
		if (mesh.getDimension() == 3)
		  gradF[2] = gradF_X[2][v];	

		const T_Scalar dFluxSecondary = gradF*secondaryCoeff;
	
		const T_Scalar dFlux = diffCoeff*((xCell[c1])[v]-(xCell[c0])[v]) + dFluxSecondary;

		(rCell[c0])[v] += dFlux;
		(rCell[c1])[v] -= dFlux;
        
		(assembler.getCoeff01(f))[v] +=diffCoeff;
		(assembler.getCoeff10(f))[v] +=diffCoeff;

		(diag[c0])[v] -= diffCoeff;
		(diag[c1])[v] -= diffCoeff;
	      }
	      
	  }
      }
  }
 private:
  const GeomFields& _geomFields;
  Field& _varField;
  const Field& _diffusivityField; 
  const Field& _varGradientField;
};

#endif
