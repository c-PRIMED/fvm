#ifndef _GRADIENTMODEL_H_
#define _GRADIENTMODEL_H_

#include "Model.h"

#include "Gradient.h"
#include "GradientMatrix.h"
#include "GeomFieldSet.h"

#include "NumType.h"
#include "GradientMatrix.h"
#include "MatrixCreator.h"

#include "Mesh.h"
#include "MeshInterface.h"

template<class X>
class GradientModel : public Model
{
public:
  typedef typename NumTypeTraits<X>::T_Scalar T_Scalar;
  
  typedef Array<X> XArray;
  typedef Gradient<X> GradType;
  typedef Array<GradType> GradArray;

  class LSQGMCalculator : public MatrixCreator
  {
  public:

    typedef GradientMatrix<T_Scalar> T_GradientMatrix;
    typedef Array<T_Scalar> TArray;
    typedef Vector<T_Scalar,3> VectorT3;
    typedef Array<VectorT3> VectorT3Array;
    typedef typename T_GradientMatrix::PairWiseAssembler GradientMatrixAssembler;

    LSQGMCalculator(const Mesh& mesh, const Field& coordField,
                    const Field& areaField, const Field& volumeField) :
      _mesh(mesh),
      _coordField(coordField),
      _areaField(areaField),
      _volumeField(volumeField)
    {}
    
    Matrix*
    getLeastSquaresGradientMatrix3D(const Mesh& mesh) const
    {
      const StorageSite& cells = mesh.getCells();
      const StorageSite& faces = mesh.getFaces();
      
      const CRConnectivity& faceCells = mesh.getConnectivity(faces,cells);
      const CRConnectivity& cellCells = mesh.getConnectivity(cells,cells);
      
      const int cellCount = cells.getSelfCount();
      const int faceCount = faces.getSelfCount();
      
      T_GradientMatrix* gMPtr =      new T_GradientMatrix(cellCells);
      T_GradientMatrix& gM = *gMPtr;
      GradientMatrixAssembler& assembler = gM.getPairWiseAssembler(faceCells);

      VectorT3Array& coeffs = gM.getCoeffs();

      const VectorT3Array& cellCentroid =
        SafeCast<VectorT3Array>(_coordField[cells]);

      const VectorT3Array& faceArea =
        SafeCast<VectorT3Array>(_areaField[faces]);

      const TArray& cellVolume = SafeCast<TArray>(_volumeField[cells]);

      const T_Scalar epsilon(1e-6);
    
      coeffs.zero();

      Array<bool> isDegenerate(cells.getCount());
      isDegenerate = false;
    
      for(int f=0; f<faceCount; f++)
      {
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          const VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];
          assembler.getCoeff01(f)=ds;
          assembler.getCoeff10(f)=-ds;
      }

      const Array<int>& row = cellCells.getRow();
    
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
                
                  coeffs[inb][0] = Kxx*ds[0] + Kxy*ds[1] + Kxz*ds[2];
                  coeffs[inb][1] = Kxy*ds[0] + Kyy*ds[1] + Kyz*ds[2];
                  coeffs[inb][2] = Kxz*ds[0] + Kyz*ds[1] + Kzz*ds[2];
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
  
    Matrix*
    getLeastSquaresGradientMatrix2D(const Mesh& mesh) const
    {
      const StorageSite& cells = mesh.getCells();
      const StorageSite& faces = mesh.getFaces();
    
      const CRConnectivity& faceCells = mesh.getConnectivity(faces,cells);
      const CRConnectivity& cellCells = mesh.getConnectivity(cells,cells);
    
      const int cellCount = cells.getSelfCount();
      const int faceCount = faces.getSelfCount();
    
      T_GradientMatrix* gMPtr =      new T_GradientMatrix(cellCells);
      T_GradientMatrix& gM = *gMPtr;
      GradientMatrixAssembler& assembler = gM.getPairWiseAssembler(faceCells);

      VectorT3Array& coeffs = gM.getCoeffs();

      const VectorT3Array& cellCentroid =
        SafeCast<VectorT3Array>(_coordField[cells]);

      const VectorT3Array& faceArea =
        SafeCast<VectorT3Array>(_areaField[faces]);

      const TArray& cellVolume = SafeCast<TArray>(_volumeField[cells]);

      const T_Scalar epsilon(1e-6);
    
      coeffs.zero();

      Array<bool> isDegenerate(cells.getCount());
      isDegenerate = false;
    
      for(int f=0; f<faceCount; f++)
      {
          const int c0 = faceCells(f,0);
          const int c1 = faceCells(f,1);
          const VectorT3 ds = cellCentroid[c1]-cellCentroid[c0];
          assembler.getCoeff01(f)=ds;
          assembler.getCoeff10(f)=-ds;
      }

      const Array<int>& row = cellCells.getRow();
    
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
                
                  coeffs[inb][0] = Kxx*ds[0] + Kxy*ds[1];
                  coeffs[inb][1] = Kxy*ds[0] + Kyy*ds[1];
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

    Matrix* create() const
    {
        const int dim = _mesh.getInt("dimension");
        if (dim==2)
          return getLeastSquaresGradientMatrix2D(_mesh);
        else
          return getLeastSquaresGradientMatrix3D(_mesh);
    }
  private:
    const Mesh& _mesh;
    const Field& _coordField;
    const Field& _areaField;
    const Field& _volumeField;
  };
  
  GradientModel(const Args& args) :
    Model(args),
    _varField(getChildRef<Field>("varField")),
    _gradientField(getChildRef<Field>("gradientField")),
    _geomFields(getChildRef<GeomFieldSet>("geomFields")),
    _coordField(_geomFields.getChildRef<Field>("coordinates")),
    _areaField(_geomFields.getChildRef<Field>("area")),
    _volumeField(_geomFields.getChildRef<Field>("volume"))
  {
    const int numMeshes = _domain.getMeshCount();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& umesh = SafeCast<Mesh>(_domain.getMesh(n));
        const StorageSite& cells = umesh.getCells();

        if (!_geomFields.hasGradientMatrixCreator(cells))
          _geomFields.addGradientMatrixCreator(cells,
                                               new LSQGMCalculator(umesh,_coordField,
                                                                   _areaField,
                                                                   _volumeField));
    }
    logCtor();
  }
                        
  virtual ~GradientModel()
  {
    const int numMeshes = _domain.getMeshCount();
    for (int n=0; n<numMeshes; n++)
    {
        const Mesh& umesh = SafeCast<Mesh>(_domain.getMesh(n));
        const StorageSite& cells = umesh.getCells();

        _geomFields.removeGradientMatrixCreator(cells);
    }
    logDtor();
  }

  DECLARE_HT("GradientModel<"+NumTypeTraits<X>::getTypeName()+">");

  PyReturn compute(PyArgsIn args)
  {
    const int nMeshes = _domain.getMeshCount();
    
    Py_BEGIN_ALLOW_THREADS
      //#pragma omp parallel for
    for(int n=0; n<nMeshes; n++)
    {
        const Mesh& mesh = SafeCast<Mesh>(_domain.getMesh(n));

        const StorageSite& cells = mesh.getCells();
        const GradientMatrix<T_Scalar>& gradMatrix =
          SafeCast<GradientMatrix<T_Scalar> > (_geomFields.getGradientMatrix(cells));
        const XArray& var = SafeCast<XArray>(_varField[cells]);
        GradArray* gradPtr = gradMatrix.getGradient(var);
        _gradientField.addArray(cells,gradPtr);


        // copy boundary values from adjacent cells
        const int nBoundaryGroups = mesh.getBoundaryGroupCount();
       
        for(int nfg=0; nfg<nBoundaryGroups; nfg++)
        {
            const Mesh::FaceGroup& fg = mesh.getBoundaryGroup(nfg);
            const StorageSite& faces = fg.site;
            const CRConnectivity& faceCells = mesh.getConnectivity(faces,cells);
            const int faceCount = faces.getCount();
            
            for(int f=0; f<faceCount; f++)
            {
                const int c0 = faceCells(f,0);
                const int c1 = faceCells(f,1);
                
                (*gradPtr)[c1] = (*gradPtr)[c0];
            }
        }

    }
    Py_END_ALLOW_THREADS
      ;

    for(int nmi=0; nmi<_domain.getMeshInterfaceCount(); nmi++)
    {
        const MeshInterface& mi = _domain.getMeshInterface(nmi);
        GradArray* gradPtr = new GradArray(mi.getCells().getCount());
        gradPtr->zero();
        _gradientField.addArray(mi.getCells(), gradPtr);
    }

    for(int n=0; n<nMeshes; n++)
    {
        const Mesh& mesh = SafeCast<Mesh>(_domain.getMesh(n));

        const StorageSite& cells = mesh.getCells();
        _gradientField.syncGather(cells);
    }
    return PyReturn();
  }

private:
  const Field& _varField;
  Field& _gradientField;
  GeomFieldSet& _geomFields;
  Field& _coordField;
  Field& _areaField;
  Field& _volumeField;
};

template<class X>
void
GradientModel<X>::addMethods()
{
  INHERIT_METHODS(Model);
  ADD_METHOD(GradientModel<X>,compute);
}

REGISTER_HT_TEMPLATE(<class X>, GradientModel, <X>);
#endif

