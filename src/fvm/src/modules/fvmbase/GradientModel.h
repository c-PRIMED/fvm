// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _GRADIENTMODEL_H_
#define _GRADIENTMODEL_H_

#include "Model.h"

#include "Gradient.h"
#include "GradientMatrix.h"
#include "GeomFields.h"

#include "NumType.h"
#include "GradientMatrix.h"
#include "SquareTensor.h"

#include "Mesh.h"


template<class T>
void
reflectGradient(Gradient<T>& gr, const Gradient<T>& g0, const Vector<T,3>& en)
{
  const T g0_dot_en_x2 = T(2.0)*(g0*en);
  for(int i=0; i<3; i++)
    gr[i] = g0[i] - g0_dot_en_x2*en[i];
}


// for a general vector treat as independent scalars
template<class T, int N>
void
reflectGradient(Gradient<Vector<T,N> >& gr,
                const Gradient<Vector<T,N> >& g0, const Vector<T,3>& en)
{
  for(int k=0; k<N; k++)
  {
      const T g0_dot_en_x2 = T(2.0)*(g0[0][k]*en[0] + g0[1][k]*en[1] + g0[2][k]*en[2]);
      for(int i=0; i<3; i++)
        gr[i][k] = g0[i][k] - g0_dot_en_x2*en[i];
  }
}

// specialization for vector 2 and 3 assume the variable is actually a
// vector in mathematical sense; this form shouldn't be used for
// independent scalars

template<class T>
void
reflectGradient(Gradient<Vector<T,2> >& gr,
                const Gradient<Vector<T,2> >& g0, const Vector<T,3>& en)
{
  const Vector<T,2> g0_dot_en_x2 = T(2.0)*(g0*en);
  for(int i=0; i<3; i++)
  {
      gr[i][0] = g0[i][0] - g0_dot_en_x2[0]*en[i];
      gr[i][1] = g0[i][1] - g0_dot_en_x2[1]*en[i];
  }
}

template<class T>
void
reflectGradient(Gradient<Vector<T,3> >& gr,
                const Gradient<Vector<T,3> >& g0, const Vector<T,3>& en)
{
  SquareTensor<T,3> R, GT0;
  
  for(int i=0;i<3;i++)
    for(int j=0; j<3; j++)
    {
        if (i==j)
          R(i,j) = 1.0 - 2*en[i]*en[j];
        else
          R(i,j) = - 2*en[i]*en[j];

        GT0(i,j) = g0[j][i];
    }

  SquareTensor<T,3> GTR(R*GT0*R);
  
  for(int i=0;i<3;i++)
    for(int j=0; j<3; j++)
    {
        gr[j][i]  = GTR(i,j);
    }
}


// non templated base class where we can store and manage the gradient matrices
// used by all templated versions
class GradientModelBase : public Model
{
public:
  GradientModelBase(const MeshList& meshes) :
    Model(meshes)
  {}

  static void clearGradientMatrix(const Mesh& mesh);
protected:
  static map<const Mesh*, shared_ptr<GradientMatrixBase> > _gradientMatricesMap;

};

template<class X>
class GradientModel : public GradientModelBase
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef Array<T_Scalar> TArray;

  typedef Array<int> IntArray;

  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<X> XArray;
  typedef Gradient<X> GradType;

  typedef Array<GradType> GradArray;

  typedef GradientMatrix<T_Scalar> GradMatrixType;
  typedef typename GradMatrixType::PairWiseAssembler GradientMatrixAssembler;

  static
  shared_ptr<GradientMatrixBase>
  getLeastSquaresGradientMatrix3D(const Mesh& mesh, const GeomFields& geomFields)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    
    const int cellCount = cells.getSelfCount();
    const int faceCount = faces.getSelfCount();
    GradMatrixType* gMPtr(new GradMatrixType(mesh));
    GradMatrixType& gM = *gMPtr;
    GradientMatrixAssembler& assembler = gM.getPairWiseAssembler(faceCells);

    const CRConnectivity& cellCells = gM.getConnectivity();
    
    VectorT3Array& coeffs = gM.getCoeffs();
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(geomFields.coordinate[cells]);

    const VectorT3Array& faceCentroid =
      dynamic_cast<const VectorT3Array&>(geomFields.coordinate[faces]);
    
    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array& >(geomFields.area[faces]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(geomFields.volume[cells]);

    const T_Scalar epsilon(1e-6);
    
    const IntArray& ibType = dynamic_cast<const IntArray&>(geomFields.ibType[cells]);
    
    coeffs.zero();

    Array<bool> isDegenerate(cells.getCount());
    isDegenerate = false;

    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];

	//fix the distance for ibfaces
	if (((ibType[c0] == Mesh::IBTYPE_FLUID)
             && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID)
             && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
	    if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
	        ds = faceCentroid[f]-cellCentroid[c0];
            }
            else
            {
	        ds = cellCentroid[c1]-faceCentroid[f];
            }
        }
        const T_Scalar dsMag = mag(ds);
        assembler.getCoeff01(f)=ds/dsMag;
        assembler.getCoeff10(f)=-ds/dsMag;
    }
  
    const Array<int>& row = cellCells.getRow();
    //const Array<int>& col = cellCells.getCol();
    
    for(int nc=0; nc<cellCount; nc++)
    {
        T_Scalar Ixx(0), Iyy(0), Izz(0);
        T_Scalar Ixy(0), Ixz(0), Iyz(0);
        
        for(int inb=row[nc]; inb<row[nc+1]; inb++)
        {
            const VectorT3& ds = coeffs[inb];
            Ixx += ds[0]*ds[0];
            Iyy += ds[1]*ds[1];
            Izz += ds[2]*ds[2];
            Ixy += ds[0]*ds[1];
            Ixz += ds[0]*ds[2];
            Iyz += ds[1]*ds[2];
        }
    

        const T_Scalar det = Ixx*(Iyy*Izz-Iyz*Iyz) -Ixy*(Ixy*Izz-Iyz*Ixz)
          + Ixz*(Ixy*Iyz-Iyy*Ixz);
        //T_Scalar det = NumTypeTraits<T_Scalar>::getZero();
        if (det > epsilon)
        {
            const T_Scalar Kxx = (Iyy*Izz-Iyz*Iyz)/det;
            const T_Scalar Kxy = -(Ixy*Izz-Iyz*Ixz)/det;
            const T_Scalar Kxz = (Ixy*Iyz-Iyy*Ixz)/det;
            const T_Scalar Kyy = (Ixx*Izz-Ixz*Ixz)/det;
            const T_Scalar Kyz = -(Ixx*Iyz-Ixy*Ixz)/det;
            const T_Scalar Kzz = (Ixx*Iyy-Ixy*Ixy)/det;
            for(int inb=row[nc]; inb<row[nc+1]; inb++)
            {
                // make a copy
                VectorT3 ds(coeffs[inb]);
                //const int j = col[inb];
		//T_Scalar dsMag = mag(cellCentroid[j]-cellCentroid[nc]);
                coeffs[inb][0] = (Kxx*ds[0] + Kxy*ds[1] + Kxz*ds[2]);
                coeffs[inb][1] = (Kxy*ds[0] + Kyy*ds[1] + Kyz*ds[2]);
                coeffs[inb][2] = (Kxz*ds[0] + Kyz*ds[1] + Kzz*ds[2]);
            }
        }
        else
        {
            isDegenerate[nc] = true;
        }
    }
    
    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];

	//fix the distance for ibfaces
	if (((ibType[c0] == Mesh::IBTYPE_FLUID)
             && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID)
             && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
	    if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
	        ds = faceCentroid[f]-cellCentroid[c0];
            }
            else
            {
	        ds = cellCentroid[c1]-faceCentroid[f];
            }
        }
        const T_Scalar dsMag = mag(ds);
        assembler.getCoeff01(f) /= dsMag;
        assembler.getCoeff10(f) /= dsMag;
    }


    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        if (isDegenerate[c0])
        {
            assembler.getCoeff01(f)= T_Scalar(0.5)*faceArea[f]/cellVolume[c0];
        }
        if (isDegenerate[c1])
        {
            assembler.getCoeff10(f)= T_Scalar(-0.5)*faceArea[f]/cellVolume[c1];
        }
    }

    return shared_ptr<GradientMatrixBase>(gMPtr);
  }

  static
  shared_ptr<GradientMatrixBase>
  getLeastSquaresGradientMatrix2D(const Mesh& mesh, const GeomFields& geomFields)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
   

    const int cellCount = cells.getSelfCount();
    const int faceCount = faces.getSelfCount();
      
    GradMatrixType* gMPtr(new GradMatrixType(mesh));
    GradMatrixType& gM = *gMPtr;
    GradientMatrixAssembler& assembler = gM.getPairWiseAssembler(faceCells);

    const CRConnectivity& cellCells = gM.getConnectivity();
    
    VectorT3Array& coeffs = gM.getCoeffs();
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(geomFields.coordinate[cells]);

    const VectorT3Array& faceCentroid =
      dynamic_cast<const VectorT3Array&>(geomFields.coordinate[faces]);

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array& >(geomFields.area[faces]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(geomFields.volume[cells]);

    const IntArray& ibType = dynamic_cast<const IntArray&>(geomFields.ibType[cells]);

    const T_Scalar epsilon(1e-26);
    
    coeffs.zero();

    Array<bool> isDegenerate(cells.getCount());
    isDegenerate = false;
    
    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];

	//fix the distance for ibfaces
	if (((ibType[c0] == Mesh::IBTYPE_FLUID)
             && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID)
             && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
	    if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
	        ds = faceCentroid[f]-cellCentroid[c0];
            }
            else
            {
	        ds = cellCentroid[c1]-faceCentroid[f];
            }
        }
        const T_Scalar dsMag = mag(ds);
        assembler.getCoeff01(f)=ds/dsMag;
        assembler.getCoeff10(f)=-ds/dsMag;
    }

    const Array<int>& row = cellCells.getRow();
    //const Array<int>& col = cellCells.getCol();
    
    for(int nc=0; nc<cellCount; nc++)
    {
        T_Scalar Ixx(0), Iyy(0);
        T_Scalar Ixy(0);
        
        for(int inb=row[nc]; inb<row[nc+1]; inb++)
        {
            const VectorT3& ds = coeffs[inb];
            Ixx += ds[0]*ds[0];
            Iyy += ds[1]*ds[1];
            Ixy += ds[0]*ds[1];
        }

        const T_Scalar det = Ixx*Iyy-Ixy*Ixy;
        //T_Scalar det = 0;
        if (det > epsilon)
        {
            const T_Scalar Kxx = Iyy/det;
            const T_Scalar Kxy = -Ixy/det;
            const T_Scalar Kyy = Ixx/det;
            for(int inb=row[nc]; inb<row[nc+1]; inb++)
            {
                // make a copy
                VectorT3 ds(coeffs[inb]);

                //const int j = col[inb];
                //T_Scalar dsMag = mag(cellCentroid[j]-cellCentroid[nc]);          
	
                coeffs[inb][0] = (Kxx*ds[0] + Kxy*ds[1]);
                coeffs[inb][1] = (Kxy*ds[0] + Kyy*ds[1]);
                coeffs[inb][2] = 0;
            }
        }
        else
        {
            isDegenerate[nc] = true;
        }
    }

    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];

	//fix the distance for ibfaces
	if (((ibType[c0] == Mesh::IBTYPE_FLUID)
             && (ibType[c1] == Mesh::IBTYPE_BOUNDARY)) ||
            ((ibType[c1] == Mesh::IBTYPE_FLUID)
             && (ibType[c0] == Mesh::IBTYPE_BOUNDARY)))
        {
	    if (ibType[c0] == Mesh::IBTYPE_FLUID)
            {
	        ds = faceCentroid[f]-cellCentroid[c0];
            }
            else
            {
	        ds = cellCentroid[c1]-faceCentroid[f];
            }
        }
        const T_Scalar dsMag = mag(ds);
        assembler.getCoeff01(f) /= dsMag;
        assembler.getCoeff10(f) /= dsMag;
    }

    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        if (isDegenerate[c0])
        {
            assembler.getCoeff01(f)= T_Scalar(0.5)*faceArea[f]/cellVolume[c0];
        }
        if (isDegenerate[c1])
        {
            assembler.getCoeff10(f)= T_Scalar(-0.5)*faceArea[f]/cellVolume[c1];
        }
    }


    return shared_ptr<GradientMatrixBase>(gMPtr);
  }

  
  GradientModel(const MeshList& meshes,
                const Field& varField, Field& gradientField,
                const GeomFields& geomFields) :
    GradientModelBase(meshes),
    _varField(varField),
    _gradientField(gradientField),
    _geomFields(geomFields)
  {
    logCtor();
  }
                        
  virtual ~GradientModel()
  {}

  static GradMatrixType&
  getGradientMatrix(const Mesh& mesh, const GeomFields& geomFields)
  {
    if (_gradientMatricesMap.find(&mesh) == _gradientMatricesMap.end())
    {
        if (mesh.getDimension() == 2)
        {
            _gradientMatricesMap[&mesh] = getLeastSquaresGradientMatrix2D(mesh,geomFields);
        }
        else
        {
            _gradientMatricesMap[&mesh] = getLeastSquaresGradientMatrix3D(mesh,geomFields);
        }
    }
    return dynamic_cast<GradMatrixType&>(*_gradientMatricesMap[&mesh]);
  }

  void init() {}
  
  void compute()
  {
 
    const int numMeshes = _meshes.size();

    for (int n=0; n<numMeshes; n++)
    {
     
        const Mesh& mesh = *_meshes[n];
	
	if (mesh.isShell() == false){
        const StorageSite& cells = mesh.getCells();
	
        GradMatrixType& gradMatrix = getGradientMatrix(mesh,_geomFields);
	
        const XArray& var = dynamic_cast<const XArray&>(_varField[cells]);
        shared_ptr<GradArray> gradPtr = gradMatrix.getGradient(var);

        _gradientField.addArray(cells,gradPtr);

        // fix values in cells adjacent to IB Faces

        const StorageSite& ibFaces = mesh.getIBFaces();
        const int nIBFaces = ibFaces.getCount();
        if (nIBFaces > 0)
        {
            const Array<int>& ibFaceList = mesh.getIBFaceList();
            
            const CRConnectivity& faceCells = mesh.getAllFaceCells();
            
            const IntArray& ibType =
              dynamic_cast<const IntArray&>(_geomFields.ibType[cells]);
            
            GradientMatrixAssembler& assembler =
              gradMatrix.getPairWiseAssembler(faceCells);
              
            const XArray& varIB =
              dynamic_cast<const XArray&>(_varField[ibFaces]);
       
            for (int nf=0; nf<nIBFaces; nf++)
            {
                const int f = ibFaceList[nf];
                const int c0 = faceCells(f,0);
                const int c1 = faceCells(f,1);
                
                if (ibType[c0]==Mesh::IBTYPE_FLUID)
                {
                    (*gradPtr)[c0].accumulate(assembler.getCoeff01(f),
                                              varIB[nf]-var[c1]);
                }
                else
                {
                    (*gradPtr)[c1].accumulate(assembler.getCoeff10(f),
                                              varIB[nf]-var[c0]);
                }
            }
        }

        // copy boundary values from adjacent cells
       
        foreach(const FaceGroupPtr fgPtr, mesh.getAllFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            const CRConnectivity& faceCells = mesh.getFaceCells(faces);
            const int faceCount = faces.getCount();
	    if ((fg.groupType!="interior") && (fg.groupType!="interface")
                && (fg.groupType!="dielectric interface"))
	    {
	        if (fg.groupType == "symmetry")
		{
		    const VectorT3Array& faceArea =
		      dynamic_cast<const VectorT3Array&>(_geomFields.area[faces]);
		    const TArray& faceAreaMag =
		      dynamic_cast<const TArray&>(_geomFields.areaMag[faces]);
		    for(int f=0; f<faceCount; f++)
		    {
		        const int c0 = faceCells(f,0);
			const int c1 = faceCells(f,1);
			const VectorT3 en = faceArea[f]/faceAreaMag[f];
			reflectGradient((*gradPtr)[c1], (*gradPtr)[c0], en);
		    }
		}
		else
		{
		    for(int f=0; f<faceCount; f++)
		    {
		        const int c0 = faceCells(f,0);
			const int c1 = faceCells(f,1);
			
			(*gradPtr)[c1] = (*gradPtr)[c0];
		    }
		}
	    }
	    }
        }
    }
    
    
    //create values
    for (int n=0; n<numMeshes; n++ ){
       const Mesh& mesh = *_meshes[n];
       if (mesh.isShell() == false){
          GradMatrixType& gradMatrix = getGradientMatrix(mesh,_geomFields);
          gradMatrix.createScatterGatherValuesBuffer();
          //partitioner interfaces
          gradMatrix.syncValues();
       }  
    }       
    
//skipping for multiple meshes, we have to fixe it,    
#if 0    
    for (int n=0; n<numMeshes; n++ ){
       const Mesh& mesh = *_meshes[n];
       if (mesh.isShell() == false){
          GradMatrixType& gradMatrix = getGradientMatrix(mesh,_geomFields);
          //mesh interfaces
          gradMatrix.recvScatterGatherValuesBufferLocal();
       }   
    }
     
#endif            
    
    typedef map<const Mesh*, shared_ptr<GradientMatrixBase> > gradType;
    foreach( const gradType::value_type& mpos, _gradientMatricesMap){
       GradientMatrixBase& gMtrxBase = *mpos.second;
       gMtrxBase.syncLocal();
    }
    _gradientField.syncLocal();

  }

private:
  const Field& _varField;
  Field& _gradientField;
  const GeomFields& _geomFields;
};


#endif

