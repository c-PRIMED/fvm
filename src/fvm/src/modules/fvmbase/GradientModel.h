#ifndef _GRADIENTMODEL_H_
#define _GRADIENTMODEL_H_

#include "Model.h"

#include "Gradient.h"
#include "GradientMatrix.h"
#include "GeomFields.h"

#include "NumType.h"
#include "GradientMatrix.h"

#include "Mesh.h"


template<class X>
class GradientModel : public Model
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef Array<T_Scalar> TArray;

  typedef Vector<T_Scalar,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  typedef Array<X> XArray;
  typedef Gradient<X> GradType;

  typedef Array<GradType> GradArray;

  typedef GradientMatrix<T_Scalar> GradMatrixType;
  typedef typename GradMatrixType::PairWiseAssembler GradientMatrixAssembler;

  static
  shared_ptr<GradMatrixType>
  getLeastSquaresGradientMatrix3D(const Mesh& mesh, const GeomFields& geomFields)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const CRConnectivity& cellCells = mesh.getCellCells();
    
    const int cellCount = cells.getSelfCount();
    const int faceCount = faces.getSelfCount();
      
    shared_ptr<GradMatrixType> gMPtr(new GradMatrixType(cellCells));
    GradMatrixType& gM = *gMPtr;
    GradientMatrixAssembler& assembler = gM.getPairWiseAssembler(faceCells);
    
    VectorT3Array& coeffs = gM.getCoeffs();
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(geomFields.coordinate[cells]);

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array& >(geomFields.area[faces]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(geomFields.volume[cells]);

    const T_Scalar epsilon(1e-6);
    
    coeffs.zero();

    Array<bool> isDegenerate(cells.getCount());
    isDegenerate = false;
    
    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        const VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];
        const T_Scalar dsMag = mag(ds);
        assembler.getCoeff01(f)=ds/dsMag;
        assembler.getCoeff10(f)=-ds/dsMag;
    }

    const Array<int>& row = cellCells.getRow();
    const Array<int>& col = cellCells.getCol();
    
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
                const int j = col[inb];
                const T_Scalar dsMag = mag(cellCentroid[j]-cellCentroid[nc]);
                
                coeffs[inb][0] = (Kxx*ds[0] + Kxy*ds[1] + Kxz*ds[2])/dsMag;
                coeffs[inb][1] = (Kxy*ds[0] + Kyy*ds[1] + Kyz*ds[2])/dsMag;
                coeffs[inb][2] = (Kxz*ds[0] + Kyz*ds[1] + Kzz*ds[2])/dsMag;
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
        if (isDegenerate[c0])
        {
            assembler.getCoeff01(f)= T_Scalar(0.5)*faceArea[f]/cellVolume[c0];
        }
        if (isDegenerate[c1])
        {
            assembler.getCoeff10(f)= T_Scalar(-0.5)*faceArea[f]/cellVolume[c1];
        }
    }

    return gMPtr;
  }

  static
  shared_ptr<GradMatrixType>
  getLeastSquaresGradientMatrix2D(const Mesh& mesh, const GeomFields& geomFields)
  {
    const StorageSite& cells = mesh.getCells();
    const StorageSite& faces = mesh.getFaces();
    
    const CRConnectivity& faceCells = mesh.getAllFaceCells();
    const CRConnectivity& cellCells = mesh.getCellCells();
    
    const int cellCount = cells.getSelfCount();
    const int faceCount = faces.getSelfCount();
      
    shared_ptr<GradMatrixType> gMPtr(new GradMatrixType(cellCells));
    GradMatrixType& gM = *gMPtr;
    GradientMatrixAssembler& assembler = gM.getPairWiseAssembler(faceCells);
    
    VectorT3Array& coeffs = gM.getCoeffs();
    
    const VectorT3Array& cellCentroid =
      dynamic_cast<const VectorT3Array&>(geomFields.coordinate[cells]);

    const VectorT3Array& faceArea =
      dynamic_cast<const VectorT3Array& >(geomFields.area[faces]);

    const TArray& cellVolume =
      dynamic_cast<const TArray&>(geomFields.volume[cells]);


    const T_Scalar epsilon(1e-26);
    
    coeffs.zero();

    Array<bool> isDegenerate(cells.getCount());
    isDegenerate = false;
    
    for(int f=0; f<faceCount; f++)
    {
        const int c0 = faceCells(f,0);
        const int c1 = faceCells(f,1);
        const VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];
        const T_Scalar dsMag = mag(ds);
        assembler.getCoeff01(f)=ds/dsMag;
        assembler.getCoeff10(f)=-ds/dsMag;
    }

    const Array<int>& row = cellCells.getRow();
    const Array<int>& col = cellCells.getCol();
    
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

                const int j = col[inb];
                const T_Scalar dsMag = mag(cellCentroid[j]-cellCentroid[nc]);
                
                coeffs[inb][0] = (Kxx*ds[0] + Kxy*ds[1])/dsMag;
                coeffs[inb][1] = (Kxy*ds[0] + Kyy*ds[1])/dsMag;
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
        if (isDegenerate[c0])
        {
            assembler.getCoeff01(f)= T_Scalar(0.5)*faceArea[f]/cellVolume[c0];
        }
        if (isDegenerate[c1])
        {
            assembler.getCoeff10(f)= T_Scalar(-0.5)*faceArea[f]/cellVolume[c1];
        }
    }

    return gMPtr;
  }

  
  GradientModel(const MeshList& meshes,
                const Field& varField, Field& gradientField,
                const GeomFields& geomFields) :
    Model(meshes),
    _varField(varField),
    _gradientField(gradientField),
    _geomFields(geomFields)
  {
    logCtor();
  }
                        
  virtual ~GradientModel()
  {}

  static const GradMatrixType&
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
    return *_gradientMatricesMap[&mesh];
  }

  void init() {}
  
  void compute()
  {
    const int numMeshes = _meshes.size();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& mesh = *_meshes[n];

        const StorageSite& cells = mesh.getCells();

        const GradMatrixType& gradMatrix = getGradientMatrix(mesh,_geomFields);

        const XArray& var = dynamic_cast<const XArray&>(_varField[cells]);
        shared_ptr<GradArray> gradPtr = gradMatrix.getGradient(var);
        _gradientField.addArray(cells,gradPtr);


        // copy boundary values from adjacent cells
       
        foreach(const FaceGroupPtr fgPtr, mesh.getBoundaryFaceGroups())
        {
            const FaceGroup& fg = *fgPtr;
            const StorageSite& faces = fg.site;
            const CRConnectivity& faceCells = mesh.getFaceCells(faces);
            const int faceCount = faces.getCount();
            
            for(int f=0; f<faceCount; f++)
            {
                const int c0 = faceCells(f,0);
                const int c1 = faceCells(f,1);
                
                (*gradPtr)[c1] = (*gradPtr)[c0];
            }
        }

    }

    _gradientField.syncLocal();
  }

private:
  const Field& _varField;
  Field& _gradientField;
  const GeomFields& _geomFields;
  static map<const Mesh*, shared_ptr<GradMatrixType> > _gradientMatricesMap;
};

template<class T>
map<const Mesh*, shared_ptr<typename GradientModel<T>::GradMatrixType> >
GradientModel<T>::_gradientMatricesMap;

#endif

