#ifndef _UMESHMETRICSCALCULATOR_H_
#define _UMESHMETRICSCALCULATOR_H_

#include "Model.h"

#include "UMesh.h"

#include "NumType.h"
#include "ArrayCreator.h"
#include "Array.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "GeomFieldSet.h"
#include "StorageSite.h"
#include "Domain.h"
#include "MeshInterface.h"

class UMeshMetricsCalculatorBase : public Model
{
public:

  static FieldLabel Coordinates;
  static FieldLabel Area;
  static FieldLabel Volume;
};

template<class T>
class UMeshMetricsCalculator : public UMeshCalculatorBase
{
public:

  virtual ArrayBase*
  calculateFaceCentroid(const Mesh& mesh) const
  {
      const StorageSite& faces = _mesh.getFaces();
      const StorageSite& nodes = _mesh.getNodes();
      const CRConnectivity& faceNodes = _mesh.getConnectivity(faces,nodes);

      const int count = faces.getCount();
      VectorT3Array *fcPtr = new VectorT3Array(count);
      VectorT3Array& fc = *fcPtr;

      const VectorT3Array& nodeCoord =
        SafeCast<VectorT3Array>(_coordField[nodes]);
      for(int f=0; f<count; f++)
      {
          const int numNodes = faceNodes.getCount(f);
          if (numNodes > 0)
          {
              fc[f] = nodeCoord[faceNodes(f,0)];
              for (int nf=1; nf<numNodes; nf++)
                fc[f] += nodeCoord[faceNodes(f,nf)];
              fc[f] /= T(numNodes);
          }
      }

      const VectorT3Array& faceArea =
        SafeCast<VectorT3Array>(_areaField[faces]);
      const TArray& faceAreaMag =
        SafeCast<TArray>(_areaMagField[faces]);
      
      // corrections for non-planar quad and polygonal faces
      const T twoThirds(2./3.);
      const T half(0.5);
      for(int f=0; f<count; f++)
      {
          const int numNodes = faceNodes.getCount(f);
          if (numNodes > 3)
          {
              const VectorT3 en = faceArea[f]/faceAreaMag[f];
              T denom(0.0);
              VectorT3 cfc(NumTypeTraits<VectorT3>::getZero());

              for (int nn=0; nn<numNodes; nn++)
              {
                  const int n0 = faceNodes(f,nn);
                  const int n1 = faceNodes(f,(nn+1)%numNodes);
                  const VectorT3 rc0 = nodeCoord[n0] - fc[f];
                  const VectorT3 rc1 = nodeCoord[n1] - fc[f];
                  const VectorT3 triArea = half*cross(rc0,rc1);
                  const T triAreaP = dot(triArea,en); 
                  VectorT3 xm = half*(nodeCoord[n0]+nodeCoord[n1]);

                  cfc += twoThirds*(xm-fc[f])*triAreaP;
                  denom += triAreaP;
              }
              cfc /= denom;
              fc[f] += cfc;
          }
      }
      
      return fcPtr;
  }
  
  virtual ArrayBase* calculateCellCentroid(const Mesh &mesh) const
  {
      const StorageSite& cells = _mesh.getCells();
      
      const int cellCount = cells.getCount();
      VectorT3Array *ccPtr = new VectorT3Array(cellCount);

      const StorageSite& faces = _mesh.getFaces();
      const StorageSite& cells = _mesh.getCells();

      const VectorT3Array& faceCentroid =
        SafeCast<VectorT3Array>(_coordField[faces]);
      const TArray& faceAreaMag = SafeCast<TArray>(_areaMagField[faces]);
      const int cellCount = cells.getCount();
      const int selfCellCount = cells.getSelfCount();
      const int faceCount = faces.getCount();
      const CRConnectivity& faceCells = _mesh.getConnectivity(faces,cells);

      {
          TArray weight(cellCount);
          
          weight.zero();
      }
      TArray weight(cellCount);

      weight.zero();
      cellCentroid.zero();
      
      for(int f=0; f<faceCount; f++)
      {
          for(int nc=0; nc<faceCells.getCount(f); nc++)
          {
              const int c = faceCells(f,nc);
              cellCentroid[c] += faceCentroid[f]*faceAreaMag[f];
              weight[c] += faceAreaMag[f];
          }
      }

      for(int c=0; c<selfCellCount; c++)
        cellCentroid[c] /= weight[c];

      const StorageSite& cells = _mesh.getCells();
      const int nBoundaryGroups = _mesh.getBoundaryGroupCount();
      for(int nfg=0; nfg<nBoundaryGroups; nfg++)
      {
          const Mesh::FaceGroup& fg = _mesh.getBoundaryGroup(nfg);
          const StorageSite& faces = fg.site;
          const VectorT3Array& faceCentroid =
            SafeCast<VectorT3Array>(_coordField[faces]);
          const CRConnectivity& faceCells = _mesh.getConnectivity(faces,cells);
          const int faceCount = faces.getCount();
          
          for(int f=0; f<faceCount; f++)
          {
              const int c1 = faceCells(f,1);
              cellCentroid[c1] = faceCentroid[f];
          }
      }
    }
    
  class FaceAreaCalculator : public ArrayCreator
  {
  public:
    FaceAreaCalculator(Field& areaField, const StorageSite& site,
                       const UMesh& mesh, const Field& coordField) :
      ArrayCreator(areaField,site),
      _mesh(mesh),
      _coordField(coordField)
    {}

    virtual ~FaceAreaCalculator() {}

    virtual ArrayBase* create() const
    {
      const StorageSite& faces = _mesh.getFaces();
      const StorageSite& nodes = _mesh.getNodes();
      const CRConnectivity& faceNodes = _mesh.getConnectivity(faces,nodes);

      const int count = faces.getCount();
      VectorT3Array *faPtr = new VectorT3Array(count);
      VectorT3Array& fa = *faPtr;

      const T half(0.5);
      const VectorT3Array& nodeCoord =
        SafeCast<VectorT3Array>(_coordField[nodes]);
      for(int f=0; f<count; f++)
      {
          const int numNodes = faceNodes.getCount(f);

          if (numNodes == 2)
          {
              const int n0 = faceNodes(f,0);
              const int n1 = faceNodes(f,1);
              VectorT3 dr = nodeCoord[n1]-nodeCoord[n0];
              fa[f][0] = dr[1];
              fa[f][1] = -dr[0];
              fa[f][2] = 0.;
          }
          else if (numNodes == 3)
          {
              const int n0 = faceNodes(f,0);
              const int n1 = faceNodes(f,1);
              const int n2 = faceNodes(f,2);
              VectorT3 dr10 = nodeCoord[n1]-nodeCoord[n0];
              VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
              fa[f] = half*cross(dr10,dr20);
          }
          else if (numNodes == 4)
          {
              const int n0 = faceNodes(f,0);
              const int n1 = faceNodes(f,1);
              const int n2 = faceNodes(f,2);
              const int n3 = faceNodes(f,3);
              VectorT3 dr20 = nodeCoord[n2]-nodeCoord[n0];
              VectorT3 dr31 = nodeCoord[n3]-nodeCoord[n1];
              fa[f] = half*cross(dr20,dr31);
          }
          else if (numNodes > 0)
          {
              fa[f].zero();
              for (int nn=0; nn<numNodes; nn++)
              {
                  const int n0 = faceNodes(f,nn);
                  const int n1 = faceNodes(f,(nn+1)%numNodes);
                  VectorT3 xm = T(0.5)*(nodeCoord[n1]+nodeCoord[n0]);
                  VectorT3 dr = (nodeCoord[n1]-nodeCoord[n0]);

                  fa[f][0] += xm[1]*dr[2];
                  fa[f][1] += xm[2]*dr[0];
                  fa[f][2] += xm[0]*dr[1];
              }
          }
      }

      
  class FaceAreaMagCalculator : public ArrayCreator
  {
  public:
    FaceAreaMagCalculator(Field& areaMagField, const StorageSite& site,
                          const UMesh& mesh, const Field& areaField) :
      ArrayCreator(areaMagField,site),
      _mesh(mesh),
      _areaField(areaField)
    {}

    virtual ~FaceAreaMagCalculator() {}

    virtual ArrayBase* create() const
    {
      const StorageSite& faces = _mesh.getFaces();
      const VectorT3Array& faceArea =
        SafeCast<VectorT3Array>(_areaField[faces]);

      const int count = faces.getCount();
      TArray *famPtr = new TArray(count);
      TArray& fam = *famPtr;

      for(int f=0; f<count; f++)
      {
          fam[f] = mag(faceArea[f]);
      }
      return famPtr;
    }
  private:
    const UMesh& _mesh;
    const Field& _areaField;
  };
      
  
  class CellVolumeCalculator : public ArrayCreator
  {
  public:
    CellVolumeCalculator(Field& volumeField, const StorageSite& site,
                         const UMesh& mesh, const Field& coordField,
                         const Field& areaField) :
      ArrayCreator(volumeField,site),
      _mesh(mesh),
      _coordField(coordField),
      _areaField(areaField)
    {}

    virtual ~CellVolumeCalculator() {}

    virtual ArrayBase* create() const
    {
      const StorageSite& faces = _mesh.getFaces();
      const StorageSite& cells = _mesh.getCells();

      const int cellCount = cells.getCount();
      TArray *vPtr = new TArray(cellCount);
      TArray& cellVolume = *vPtr;
      
      const VectorT3Array& faceCentroid = SafeCast<VectorT3Array>(_coordField[faces]);
      const VectorT3Array& cellCentroid = SafeCast<VectorT3Array>(_coordField[cells]);
      const VectorT3Array& faceArea = SafeCast<VectorT3Array>(_areaField[faces]);
      const CRConnectivity& faceCells = _mesh.getConnectivity(faces,cells);
      
      const T dim(_mesh.getInt("dimension"));
      const int faceCount = faces.getCount();

      cellVolume.zero();
      
      for(int f=0; f<faceCount; f++)
      {
          const int c0 = faceCells(f,0);
          cellVolume[c0] += dot(faceCentroid[f]-cellCentroid[c0],faceArea[f])/dim;
          const int c1 = faceCells(f,1);
          cellVolume[c1] -= dot(faceCentroid[f]-cellCentroid[c1],faceArea[f])/dim;
      }
    

      // update boundary cells with adjacent interior cells values
      const int nBoundaryGroups = _mesh.getBoundaryGroupCount();
      for(int nfg=0; nfg<nBoundaryGroups; nfg++)
      {
          const Mesh::FaceGroup& fg = _mesh.getBoundaryGroup(nfg);
          const StorageSite& faces = fg.site;
          const CRConnectivity& faceCells = _mesh.getConnectivity(faces,cells);
          const int faceCount = faces.getCount();
          
          for(int f=0; f<faceCount; f++)
          {
              const int c0 = faceCells(f,0);
              const int c1 = faceCells(f,1);
              cellVolume[c1] = cellVolume[c0];
          }
      }
      return vPtr;
    }
  private:
    const UMesh& _mesh;
    const Field& _coordField;
    const Field& _areaField;
  };
    
  class InterfaceVolumeCalculator : public ArrayCreator
  {
  public:
    InterfaceVolumeCalculator(Field& volumeField, const StorageSite& site,
                              const MeshInterface& mi) :
      ArrayCreator(volumeField,site),
      _mi(mi)
    {}

    virtual ~InterfaceVolumeCalculator() {}

    virtual ArrayBase* create() const
    {
      const int count = _mi.getCells().getCount();
      TArray *volPtr = new TArray(count);
      volPtr->zero();
      return volPtr;
    }
    
  private:
    const MeshInterface& _mi;
  };
      
  UMeshMetricsCalculator(const Args& args) :
    Model(args),
    _geomFields(getChildRef<GeomFieldSet>("geomFields")),
    _coordField(_geomFields.getChildRef<Field>("coordinates")),
    _areaField(_geomFields.getChildRef<Field>("area")),
    _areaMagField(_geomFields.getChildRef<Field>("areaMagnitude")),
    _volumeField(_geomFields.getChildRef<Field>("volume"))
  {
    const int numMeshes = _domain.getMeshCount();
    for (int n=0; n<numMeshes; n++)
    {
        const UMesh& umesh = SafeCast<UMesh>(_domain.getMesh(n));
        const StorageSite& faces = umesh.getFaces();
        const StorageSite& cells = umesh.getCells();
        addCreator(new FaceCentroidCalculator(_coordField,faces, umesh,
                                              _areaField,_areaMagField));
        addCreator(new FaceAreaCalculator(_areaField,faces,umesh,_coordField));
        addCreator(new FaceAreaMagCalculator(_areaMagField, faces,umesh,_areaField));
        
        addCreator(new CellCentroidCalculator(_coordField,cells, umesh,
                                              _areaMagField));
        addCreator(new CellVolumeCalculator(_volumeField,cells,umesh,_coordField,
                                            _areaField));

    }

    for(int nmi=0; nmi<_domain.getMeshInterfaceCount(); nmi++)
    {
        const MeshInterface& mi = _domain.getMeshInterface(nmi);
        addCreator(new InterfaceCentroidCalculator(_coordField, mi.getCells(),mi));
        addCreator(new InterfaceVolumeCalculator(_volumeField,mi.getCells(),mi));
    }
    logCtor();
  }

  virtual ~UMeshMetricsCalculator()
  {
    logDtor();
  }

  DECLARE_HT("UMeshMetricsCalculator<"+NumTypeTraits<T>::getTypeName()+">");

  
private:
  GeomFieldSet& _geomFields;
  Field& _coordField;
  Field& _areaField;
  Field& _areaMagField;
  Field& _volumeField;
};

template<class T>
void
UMeshMetricsCalculator<T>::addMethods()
{
  INHERIT_METHODS(Model);
}


REGISTER_HT_TEMPLATE(<class T>, UMeshMetricsCalculator, <T>);
#endif
