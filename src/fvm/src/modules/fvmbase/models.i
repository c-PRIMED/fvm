
%include "std_string.i"
%include "std_vector.i"
%include "std_except.i"
%include "std_map.i"
%include "shared_ptr.i"

using namespace std;

%import "ArrayBase.i"
%import "Field.i"
%import "Mesh.i"
%import "LinearSolver.i"
%import "Vector.i"


%include "GeomFields.h"
%include "FlowFields.h"
%include "ThermalFields.h"
%include "ElectricFields.h"
%include "StructureFields.h"
%include "atype.i"
%include "Model.i"

%include "MeshMetricsCalculator.i"

%include "ThermalModel.i"

%include "ElectricModel.i"

%include "FlowModel.i"

%include "IdealGasDensityModel.i"
%include "RosselandModel.i"
%include "StructureModel.i"

 //%include "CartesianElasticity.i"

%include "MovingMeshModel.i"


%include "OneDConduction.i"

 //%include "ShockTube.i"

%include "StructureDeformationModel.i"


#ifdef USING_ATYPE_TANGENT

typedef Vector<Tangent,3> VecTangent3;
%template(VecTangent3) Vector<Tangent,3>;

#endif



%inline{
#ifdef USING_ATYPE_PC

typedef Vector< ATYPE, 3 > VecPC3;
%template(VecPC3) Vector< ATYPE, 3 >;
    boost::shared_ptr<ArrayBase> getStdDev(boost::shared_ptr<ArrayBase> abase)
    {
      Array< ATYPE >& a = dynamic_cast<Array< ATYPE >& >(*abase);
      const int length = a.getLength();
      Array<double> *sptr = new Array<double>(length);
      for(int i=0; i<length; i++)
      {
          (*sptr)[i] = a[i].stdDev();
      }
      return  boost::shared_ptr<ArrayBase > (sptr);

    }
    boost::shared_ptr<ArrayBase> getStdDevVec3(boost::shared_ptr<ArrayBase> abase)
    {
        Array< Vector<ATYPE,3> >& a = dynamic_cast<Array< Vector<ATYPE,3> >& >(*abase);
        const int length = a.getLength();
        Array<Vector<double,3> > *sptr = new Array<Vector<double,3> >(length);
        for(int i=0; i<length; i++)
          for(int j=0; j<3; j++)
          {
              (*sptr)[i][j] = a[i][j].stdDev();
          }
        return  boost::shared_ptr<ArrayBase > (sptr);
        
    }
#endif

    boost::shared_ptr<ArrayBase> newVec3Array(const int size)
    {
      return shared_ptr<ArrayBase>(new Array< Vector<ATYPE,3> >(size));
    }

}

