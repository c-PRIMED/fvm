
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
%include "VacancyFields.h"
%include "FractureFields.h"
%include "SpeciesFields.h"
%include "BatteryFields.h"



%include "atype.i"

%include "MeshMetricsCalculator.i"
%include "ThermalModel.i"
%include "VacancyModel.i"
%include "FractureModel.i"
%include "SpeciesModel.i"
%include "BatteryModel.i"




#ifndef USING_ATYPE_PC

%include "ElectricFields.h"
%include "StructureFields.h"
%include "ElectricModel.i"
%include "MovingMeshModel.i"
%include "StructureDeformationModel.i"
%include "StructureModel.i"
%include "PlateFields.h"
%include "PlateModel.i"
%include "PlateDeformationModel.i"
#endif

%include "FlowModel.i"

%include "IdealGasDensityModel.i"


 //%include "OneDConduction.i"

 //%include "ShockTube.i"



#ifdef USING_ATYPE_TANGENT

typedef Vector<Tangent,3> VecTangent3;
%template(VecTangent3) Vector<Tangent,3>;

#endif

#ifdef USING_ATYPE_PC

typedef Vector< ATYPE, 3 > VecPC3;
%template(VecPC3) Vector< ATYPE, 3 >;
#endif

%inline{

#ifdef USING_ATYPE_PC

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


