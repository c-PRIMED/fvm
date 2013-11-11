// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _BATTERYFIXINTERFACEGHOST_H_
#define _BATTERYFIXINTERFACEGHOST_H_ 

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




template<class X, class Diag, class OffDiag>
  class BatteryFixInterfaceGhost
{
 public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  typedef CRMatrix<Diag,OffDiag,X> CCMatrix;
  typedef typename CCMatrix::DiagArray DiagArray;
  typedef Array<T_Scalar> TArray;
  typedef Array<X> XArray;
  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;
 

 BatteryFixInterfaceGhost(const MeshList& meshes,
                          GeomFields& geomFields,
                          Field& varField,
                          Field& diffusivityField):
    _meshes(meshes),
    _geomFields(geomFields),
    _varField(varField),
    _diffusivityField(diffusivityField)
    {}


  void fixInterfaces()
  {

    // this function edits the interface ghost cells when an unconnected 
    // doubleShell is present at that interface.  It does 3 things:
    //
    // 1) changes the interface ghost cell's centroid to that of the interface 
    //    instead of the interior cell of the adjacent mesh to which it is connected
    //    through the scatter/gather map that we want to remain intact.
    // 2) copies the diffusivity from the interior cell to the interface ghost cell 
    //    instead of the diffusivity that the ghost cell picked up from the adjacent mesh
    // 3) copies in interface values being stored in the doubleShell mesh to the 
    //    interface ghost cells instead of them holding the value from the interior 
    //    cell of the adjacent mesh.

    cout << "BATTERY FIX INTEFACE GHOST" << endl;

    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
      {
	const Mesh& mesh = *_meshes[n];	

	//only do process if it is an unconnected doubleShell
	if ((mesh.isDoubleShell())&&(!(mesh.isConnectedShell())))
	  {
	    const StorageSite& cells = mesh.getCells();
	    const CRConnectivity& cellCells = mesh.getCellCells();
	    const XArray& varCell = dynamic_cast<const XArray&>(_varField[cells]);

	    const int parentMeshID = mesh.getParentMeshID();
	    const Mesh& parentMesh = *_meshes[parentMeshID];
	    const StorageSite& parentFaces = mesh.getParentFaceGroupSite();
	    const CRConnectivity& parentFaceCells = parentMesh.getFaceCells(parentFaces);
	    const StorageSite& parentCells = parentMesh.getCells();
	    VectorT3Array& parentCellCentroid = dynamic_cast<VectorT3Array&>(_geomFields.coordinate[parentCells]);
	    VectorT3Array& parentFaceCentroid = dynamic_cast<VectorT3Array&>(_geomFields.coordinate[parentFaces]);
	    TArray& diffusivityParent = dynamic_cast<TArray&>(_diffusivityField[parentCells]);
	    XArray& varCellParent = dynamic_cast<XArray&>(_varField[parentCells]);


	    const int otherMeshID = mesh.getOtherMeshID();
	    const Mesh& otherMesh = *_meshes[otherMeshID];
	    const StorageSite& otherFaces = mesh.getOtherFaceGroupSite();
	    const CRConnectivity& otherFaceCells = otherMesh.getFaceCells(otherFaces);
	    const StorageSite& otherCells = otherMesh.getCells();
	    VectorT3Array& otherCellCentroid = dynamic_cast<VectorT3Array&>(_geomFields.coordinate[otherCells]);
	    VectorT3Array& otherFaceCentroid = dynamic_cast<VectorT3Array&>(_geomFields.coordinate[otherFaces]);
	    TArray& diffusivityOther = dynamic_cast<TArray&>(_diffusivityField[otherCells]);
	    XArray& varCellOther = dynamic_cast<XArray&>(_varField[otherCells]);


	    for (int f=0; f<parentFaces.getCount(); f++)
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
		//cout << "Parent: " << c1p << " " << (parentCellCentroid[c1p])[0] << " " << (parentFaceCentroid[f])[0] << endl;

		// change centroid
		parentCellCentroid[c1p] = parentFaceCentroid[f];

		// copy diffusivity from interior cell to ghost 
		// (instead of diff picked up from other mesh during sync)
		diffusivityParent[c1p] = diffusivityParent[c0p];

		// copy solution variable value from shell to ghost cell
		const int c0 = f;
		varCellParent[c1p] = varCell[c0];

	      }

	    for (int f=0; f<otherFaces.getCount(); f++)
	      {
		//get other mesh fluxes and coeffs
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

		//cout << "Other: " << c1o << " " << (otherCellCentroid[c1o])[0] << " " << (otherFaceCentroid[f])[0] << endl;

		// change centroid
		otherCellCentroid[c1o] = otherFaceCentroid[f];

		// copy diffusivity from interior cell to ghost 
		// (instead of diff picked up from other mesh during sync)
		diffusivityOther[c1o] = diffusivityOther[c0o];

		// copy solution variable value from shell to ghost cell
		const int c1 = cellCells(f,0);
		varCellOther[c1o] = varCell[c1];
	      }

	  }
      }

    //output for two material 54 cell case to check centroid changes
    if (0)
      {
	for (int n=0; n<numMeshes; n++)
	  {
	    const Mesh& mesh = *_meshes[n];	
	    const StorageSite& cells = mesh.getCells();

	    const int nCells = cells.getCount();	

	    VectorT3Array& cellCentroid = dynamic_cast<VectorT3Array&>(_geomFields.coordinate[cells]);

	    cout << "Mesh: " << n << endl;
	    for (int c=0; c<nCells; c++)
	      {
		if (((cellCentroid[c])[2] < -3.33)&&((cellCentroid[c])[2] > -3.34))
		  cout << c << ": " << (cellCentroid[c])[0] << " " << (cellCentroid[c])[1] << endl;
	      }
	  }
      }

  }
 private:
  
  GeomFields& _geomFields;
  Field& _varField;
  Field& _diffusivityField;
  const MeshList& _meshes;
  
};

#endif
