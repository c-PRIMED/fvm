#ifndef _THERMALMODEL_H_
#define _THERMALMODEL_H_

#include "Model.h"

#include "UMesh.h"

#include "NumType.h"
#include "ArrayCreator.h"
#include "Array.h"
#include "Field.h"
#include "CRConnectivity.h"
#include "FieldSet.h"
#include "StorageSite.h"
#include "Domain.h"
#include "FieldArrayInitializer.h"
#include "ZeroArrayInitializer.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "FluxJacobianMatrix.h"
#include "MeshInterface.h"
#include "Topology.h"
#include "DiagonalMatrix.h"
#include "PropertySet.h"
#include "GenericBCS.h"
#include "Vector.h"

struct ThermalFields
{
  static FieldLabel temperature;
  static FieldLabel heatFlux;
  static FieldLabel enthalpy;
  
};

template<class T>
class ThermalModel : public Model
{
public:

  typedef Array<T> TArray;
  typedef Vector<T,3> VectorT3;
  typedef Array<VectorT3> VectorT3Array;

  ThermalModel(const MeshList& meshes) :
    Model(meshes),
    _areaMagField(_geomFields.getChildRef<Field>("areaMagnitude")),
    _thermalFields(getChildRef<FieldSet>("thermalFields")),
    _temperature(createGlobalField(ThermalFields::temperature))
    _temperatureGradient(_thermalFields.getChildRef<Field>("temperatureGradient")),
    _heatFlux(_thermalFields.getChildRef<Field>("heatFlux"))
  {
    const int numMeshes = _domain.getMeshCount();

    for (int n=0; n<numMeshes; n++)
    {
        const UMesh& umesh = SafeCast<UMesh>(_domain.getMesh(n));
        const StorageSite& cells = umesh.getCells();
        addCreator(new FieldArrayInitializer<T>(_temperature, cells, *this,
                                                "initialTemperature"));

        for(int nfg=0; nfg<umesh.getFaceGroupCount(); nfg++)
        {
            const Mesh::FaceGroup& fg = umesh.getFaceGroup(nfg);
            const string groupType = fg.getString("groupType");
            if (groupType != "interior")
            {
                const StorageSite& fgSite = fg.getChildRef<StorageSite>("sSite");
                addCreator(new ZeroArrayInitializer<T>(_heatFlux,fgSite));
            }
        }
    }

    for(int nmi=0; nmi<_domain.getMeshInterfaceCount(); nmi++)
    {
        const MeshInterface& mi = _domain.getMeshInterface(nmi);
        addCreator(new FieldArrayInitializer<T>(_temperature, mi.getCells(),
                                                *this,"initialTemperature"));
    }
    logCtor();
  }

  
  virtual ~ThermalModel()
  {
    logDtor();
  }

  DECLARE_HT("ThermalModel<"+NumTypeTraits<T>::getTypeName()+">");

  PyReturn linearizeBCs(PyArgsIn args)
  {
    MultiFieldMatrix& matrix = args.getRef<MultiFieldMatrix>(0);
    MultiField& x = args.getRef<MultiField>(1);
    MultiField& b = args.getRef<MultiField>(2);

    const int nMeshes = _domain.getMeshCount();
    
    Py_BEGIN_ALLOW_THREADS
      //#pragma omp parallel for
    for(int n=0; n<nMeshes; n++)
    {
        const UMesh& mesh = SafeCast<UMesh>(_domain.getMesh(n));

        for(int nfg=0; nfg<mesh.getBoundaryGroupCount(); nfg++)
        {
            const Mesh::FaceGroup& fg = mesh.getBoundaryGroup(nfg);
            const Surface& surface = _domain.getSurface(fg);
            const StorageSite& faces = fg.getChildRef<StorageSite>("sSite");

            const BCSet& bcSet = surface.getBCSet(_thermalFields);
            const PropertySet& bcProps = bcSet.getProperties();

            const string bcType = bcSet.getBCType();

            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,_temperature,_heatFlux,
                                  matrix, x, b);

            if (bcType == "SpecifiedTemperatureWall")
            {
                const T bT(bcProps.getVar("wallTemperature"));
                gbc.applyDirichletBC(bT);
            }
            else if (bcType == "AdiabaticWall")
            {
                T zero(NumTypeTraits<T>::getZero());
                gbc.applyNeumannBC(zero);
            }
            else if (bcType == "SpecifiedHeatFluxWall")
            {
                const T specifiedFlux(bcProps.getVar("wallHeatFlux"));
                gbc.applyNeumannBC(specifiedFlux);
            }
            else if ((bcType == "ExternalConvectionWall"))
            {
                const T h(bcProps.getVar("externalHTC"));
                const T Tinf(bcProps.getVar("externalTemperature"));
                gbc.applyConvectionBC(h,Tinf);
            }
            else if ((bcType == "CoupledWall"))
            {
                gbc.applyInterfaceBC();
            }
            else if ((bcType == "FlowBoundary"))
            {
                if (!isVarNone("flowFields"))
                {
                    const FieldSet& flowFields = getChildRef<FieldSet>("flowFields");
                    const Field& convectingFluxField = flowFields.getChildRef<Field>("massFlux");
                    if (convectingFluxField.hasArray(faces))
                    {
                        const T bT(bcProps.getVar("boundaryTemperature"));
                        
                        const TArray& convFlux = SafeCast<TArray>(convectingFluxField[faces]);
                        gbc.applyFlowBC(convFlux,bT);
                    }
                    else
                      throw CException(bcType + " only valid at flow boundaries");
                }
                else
                  throw CException(bcType + " only valid if flow model is defined");
            }
            else
              throw CException(bcType + " not implemented for ThermalModel");
        }

        // interfaces
        for(int nfg=0; nfg<mesh.getInterfaceGroupCount(); nfg++)
        {
            const Mesh::FaceGroup& fg = mesh.getInterfaceGroup(nfg);
            const StorageSite& faces = fg.getChildRef<StorageSite>("sSite");
            GenericBCS<T,T,T> gbc(faces,mesh,
                                  _geomFields,_temperature,_heatFlux,
                                  matrix, x, b);
            gbc.applyInterfaceBC();
        }

    }
    Py_END_ALLOW_THREADS
      ;
      
    for(int nmi=0; nmi<_domain.getMeshInterfaceCount(); nmi++)
    {
        const MeshInterface& mi = _domain.getMeshInterface(nmi);
        applyCoupling<T,T,T>(mi,matrix, _temperature, _heatFlux, x, b);
    }
    return PyReturn();
  }
private:
  const FieldSet& _geomFields;
  const Field& _areaMagField;

  const FieldSet& _thermalFields;
  
  Field& _temperature;
  Field& _temperatureGradient;
  Field& _heatFlux;
};

template<class T>
void
ThermalModel<T>::addMethods()
{
  INHERIT_METHODS(Model);
  ADD_METHOD(ThermalModel<T>,linearizeBCs);
}

REGISTER_HT_TEMPLATE(<class T>, ThermalModel, <T>);
#endif

