// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _COMETBOUNDARY_H_
#define _COMETBOUNDARY_H_

#include "Mesh.h"
#include <math.h>
#include "NumType.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "StorageSite.h"
#include "CRConnectivity.h"
#include "GeomFields.h"
#include "pmode.h"
#include "Kspace.h"
#include "kvol.h"
#include "SquareTensor.h"
#include "GradientModel.h"
#include "FluxLimiters.h"

template<class T>
class COMETBoundary
{
 public :
  
  typedef typename NumTypeTraits<T>::T_Scalar T_Scalar;
  typedef Array<int> IntArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;  
  typedef Array<VectorT3> VectorT3Array;
  typedef Kspace<T> Tkspace;
  typedef kvol<T> Tkvol;
  typedef pmode<T> Tmode;
  typedef typename Tmode::Refl_pair Refl_pair;
  typedef Gradient<T> GradType;
  typedef Array<GradType> GradArray;
  typedef GradientModel<T> GradModelType;
  typedef typename GradModelType::GradMatrixType GradMatrix;

 COMETBoundary(const StorageSite& faces,
		const Mesh& mesh,
		const GeomFields& geomFields,
		Tkspace& kspace,
		COMETModelOptions<T>& opts,
		const int fg_id):
  _faces(faces),
    _cells(mesh.getCells()),
    _faceCells(mesh.getFaceCells(_faces)),
    _cellCells(mesh.getCellCells()),
    _faceArea(dynamic_cast<const VectorT3Array&>(geomFields.area[_faces])),
    _cellCoords(dynamic_cast<const VectorT3Array&>(geomFields.coordinate[_cells])),
    _kspace(kspace), 
    _eArray(kspace.geteArray()),
    _geomFields(geomFields),
    _mesh(mesh)
      {}

  void applyTemperatureWallFine(FloatValEvaluator<T>& bTemp) const
  {

    const GradMatrix& gradMatrix=GradModelType::getGradientMatrix(_mesh,_geomFields);

    for (int j=0; j<_faces.getCount();j++)
      {
	applyTemperatureWallFine(j, bTemp[j], gradMatrix);
      }
  }

  void applyTemperatureWallCoarse(FloatValEvaluator<T>& bTemp) const
  {
    for (int j=0; j<_faces.getCount();j++)
      {
	applyTemperatureWallCoarse(j, bTemp[j]);
      }
  }

  void applyTemperatureWallFine(int f, const T Twall, const GradMatrix& gMat) const
  {
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); 

    const int neibcount=_cellCells.getCount(c0);

    const VectorT3Array& faceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[_faces]);

    GradArray Grads(_kspace.gettotmodes());
    Grads.zero();
    TArray pointMin(_kspace.gettotmodes());
    pointMin=-1;
    TArray pointMax(_kspace.gettotmodes());
    pointMax.zero();
    TArray pointLim(_kspace.gettotmodes());
    pointLim=100;

    VectorT3 Gcoeff;

    for(int j=0;j<neibcount;j++)    //first loop to get grad and max/min vals
      {
	const int cell1=_cellCells(c0,j);
	const VectorT3 Gcoeff=gMat.getCoeff(c0, cell1);

	int c0ind=_kspace.getGlobalIndex(c0,0);
	int c1ind=_kspace.getGlobalIndex(cell1,0);
	
	for(int k=0;k<_kspace.gettotmodes();k++)
	  {
	    const T e1=_eArray[c1ind];
	    const T e0=_eArray[c0ind];
	    T& curMax=pointMax[k];
	    T& curMin=pointMin[k];
	    Grads[k].accumulate(Gcoeff, e1-e0);

	    if(e1>curMax)
	      curMax=e1;
	    
	    if(curMin==-1)
	      curMin=e1;
	    else
	      {
		if(e1<curMin)
		  curMin=e1;
	      }

	    c0ind++;
	    c1ind++;
	  }  
      }

    const StorageSite& allFaces(_mesh.getFaces());
    const VectorT3Array& allFaceCoords=
      dynamic_cast<const VectorT3Array&>(_geomFields.coordinate[allFaces]);
    const CRConnectivity& cellFaces(_mesh.getCellFaces());

    for(int j=0;j<neibcount;j++)    //second loop to calculate limiting coeff
      {
	const int f1(cellFaces(c0,j));
	const VectorT3 dr0(allFaceCoords[f1]-_cellCoords[c0]);
	int c0ind=_kspace.getGlobalIndex(c0,0);
	//vanLeer lf;
	
	for(int k=0;k<_kspace.gettotmodes();k++)
	  {
	    const T minVal=pointMin[k];
	    const T maxVal=pointMax[k];
	    const T de0(Grads[k]*dr0);
	    T& cl=pointLim[k];
	    //computeLimitCoeff(cl, _eArray[c0ind], de0, minVal, maxVal, lf);
	    c0ind++;
	  }  
      }
    
    int numK=_kspace.getlength();

    VectorT3 rVec=_cellCoords[c1]-_cellCoords[c0];
    VectorT3 fVec=faceCoords[f]-_cellCoords[c0];

    for (int k=0;k<numK;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	int numM=kv.getmodenum();
	
	for (int m=0;m<numM;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    VectorT3 vg = mode.getv();
	    const int index=mode.getIndex()-1;
	    const VectorT3 en = _faceArea[f];
	    const T vg_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
	    const int c0ind=_kspace.getGlobalIndex(c0,index);
	    const int c1ind=_kspace.getGlobalIndex(c1,index);
	    const GradType& grad=Grads[index];
	    
	    if (vg_dot_en > T_Scalar(0.0))
	      {
		const T SOU=(fVec[0]*grad[0]+fVec[1]*grad[1]+
			   fVec[2]*grad[2])*pointLim[index];
		_eArray[c1ind]=_eArray[c0ind]+SOU;
	      }
	    else
	      _eArray[c1ind]=mode.calce0(Twall);
	  }	
      }
  }

  void applyTemperatureWallCoarse(int f,const T Twall) const
  {
    
    const int c0 = _faceCells(f,0);
    const int c1 = _faceCells(f,1); 
    
    int numK=_kspace.getlength();
    
    for (int k=0;k<numK;k++)
      {
	Tkvol& kv=_kspace.getkvol(k);
	int numM=kv.getmodenum();
	
	for (int m=0;m<numM;m++)
	  {
	    Tmode& mode=kv.getmode(m);
	    VectorT3 vg = mode.getv();
	    const int index=mode.getIndex()-1;
	    const VectorT3 en = _faceArea[f];
	    const T vg_dot_en = vg[0]*en[0]+vg[1]*en[1]+vg[2]*en[2];
	    const int c0ind=_kspace.getGlobalIndex(c0,index);
	    const int c1ind=_kspace.getGlobalIndex(c1,index);
	    
	    if (vg_dot_en > T_Scalar(0.0))
	      {
		_eArray[c1ind]=_eArray[c0ind];
	      }
	    else
	      _eArray[c1ind]=mode.calce0(Twall);
	  }	
      }
  }
  
 protected:
  
  const StorageSite& _faces;
  const StorageSite& _cells;
  const CRConnectivity& _faceCells;
  const CRConnectivity& _cellCells;
  const VectorT3Array& _faceArea;
  const VectorT3Array& _cellCoords;
  Tkspace& _kspace;
  TArray& _eArray;
  const GeomFields& _geomFields;
  const Mesh& _mesh;
};

#endif
