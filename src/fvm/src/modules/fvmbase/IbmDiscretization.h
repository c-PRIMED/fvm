// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _IBMDISCRETIZATION_H_
#define _IBMDISCRETIZATION_H_

#include "CRMatrix.h"
#include "Field.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "Mesh.h"
#include "Discretization.h"
#include "StorageSite.h"
#include "DiagonalMatrix.h"
#include "DiagonalTensor.h"
#include "FlowFields.h"
#include "GeomFields.h"
#include "Vector.h"
#include "NumType.h"

#include <iostream>
#include <fstream>


  
template <class X, class Diag, class OffDiag>
class IbmDiscretization: public Discretization
{
public:

  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef typename CCMatrix::OffDiagArray OffDiagArray;

  typedef Array<X> XArray;
  typedef Array<T_Scalar> TArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  IbmDiscretization (const MeshList& meshes,
		     const GeomFields& geomFields,
		     FlowFields& flowFields) :
    Discretization(meshes),
    _geomFields(geomFields),
    _flowFields(flowFields)
{}

  void discretize(const Mesh& mesh, MultiFieldMatrix& mfmatrix,
                  MultiField& xField, MultiField& rField)
  {
    const StorageSite& cells = mesh.getCells();

    const MultiField::ArrayIndex vIndex(&_flowFields.velocity,&cells);
    CCMatrix& matrix = dynamic_cast<CCMatrix&>(mfmatrix.getMatrix(vIndex,vIndex));
    DiagArray& diag = matrix.getDiag();
    OffDiagArray& offdiag = matrix.getOffDiag();
    
    const CRConnectivity& conn=matrix.getConnectivity();
    const Array<int>& row=conn.getRow();
    const Array<int>& col=conn.getCol();
    const VectorT3Array& cellCentroid = dynamic_cast<const VectorT3Array&> (_geomFields.coordinate[cells]);
    VectorT3Array& cellVelocity = dynamic_cast<VectorT3Array&> (_flowFields.velocity[cells]);
    const XArray& xCell = dynamic_cast<const XArray&>(xField[vIndex]);
    XArray& rCell = dynamic_cast<XArray&>(rField[vIndex]); 
    
    const int nCells = cells.getSelfCount();
    const VectorT3 zeroVelocity(NumTypeTraits<VectorT3>::getZero());
    VectorT3 dR;        //local variable dR=r[c]-r[j]
    //set up solid cell veloicty to be solidVelocity
    //ofstream fp1, fp2;
    //fp1.open ("/home/linsun/Work/prism/app-memosa/src/fvm/test/ibm_velocity");
    //fp2.open ("/home/linsun/Work/prism/app-memosa/src/fvm/test/solid_velocity");
    
    // fp=fopen("/home/linsun/Work/prism/app-memosa/src/fvm/test/ibm_velocity","w");
    for (int c=0; c<nCells; c++)
    {
      ////////need to define -1 and 0 to make sense for tangent type
      int cellIBType = mesh.getIBTypeForCell(c);
      if (cellIBType == 2)               //solid cells
      {
	//cout << "cellID" << c << "velocity" << _flowFields.velocity << endl; 
	//printf("cellID = %i and velocity is %f\t%f\t%f\n", c,cellVelocity[c][0],cellVelocity[c][1],cellVelocity[c][2]); 
	diag[c] = (-1);   //ap=(-1)
	//diag[c] = NumTypeTraits<Diag>::getNegativeUnity;
	cellVelocity[c]=zeroVelocity;     //Cell velocity is set as prescribed
	for (int nbpos=row[c]; nbpos<row[c+1]; nbpos++){ 
	  //const int j=col[nbpos];
	  offdiag[nbpos]=0;              //anb=0
	  //offdiag[nbpos] = NumTypeTraits<offDiag>::getZero;
	}
	rCell[c]=0;                     //residual is zero
	//rCell[c]= NumTypeTraits<VectorT3>::getZero;
	//fprintf(fp, "%i\t%f\t%f\t%f\n", c,cellVelocity[c][0],cellVelocity[c][1],cellVelocity[c][2]);
	//fp2 << c << cellVelocity[c] << endl;
 
      }
      //fclose(fp);
     
      
      #if 1
      else if (cellIBType == 1)          //immersed boundary cells
      { 
	//cout << "cellID" << c << "velocity" << _flowFields.velocity\n; 
	diag[c] = (-1);   //aP=-1
	//diag[c] = NumTypeTraits<Diag>::getNegativeUnity;
	double sumCellDist=0;
	//offdiag[nb]=alfa[nb] (distance weight function)
	//r=-v[c]+sum(alfa[nb]*v[nb])
	/*
	for (int nbpos=row[c]; nbpos<row[c+1]; nbpos++){
	  const int j=col[nbpos];
	  dR=cellCentroid[c]-cellCentroid[j];
	  offdiag[nbpos]=1./mag(dR);
	  sumCellDist+=offdiag[nbpos];
	}
	*/
	
	for (int nbpos=row[c]; nbpos<row[c+1]; nbpos++){
	  //offdiag[nbpos]/=sumCellDist;
	  const int j=col[nbpos];
	  offdiag[nbpos]=1.0/(row[c+1]-row[c]);
	  rCell[c]-=offdiag[nbpos]*cellVelocity[j];
	}
	rCell[c]+=cellVelocity[c];
	//fp1 << c << cellVelocity[c] << endl;

      }
      #endif
    }
 //fp1.close();
 //fp2.close();
  }

  private:
  FlowFields& _flowFields;
  const GeomFields& _geomFields;
};

#endif
