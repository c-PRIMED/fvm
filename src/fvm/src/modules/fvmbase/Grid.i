%{
#include "Grid.h"
#include "Array.h"
#include "StorageSite.h"
#include "Vector.h"
#include "Field.h"
#include "CRConnectivity.h"

  %}

using namespace std;

class Grid
{
public: 
  Grid( GeomFields& geomFields, FlowFields& flowFields,
        string coordFileName, string velocityFileName);
  
  ~Grid();
  

  const StorageSite& getNodes();
  const StorageSite& getCells();
  
  void setConnFaceToGrid(Mesh& mesh, const StorageSite& faces );
  
  boost::shared_ptr<ArrayBase>
  computeInterpolatedVelocity(const StorageSite& faces);

};



 

