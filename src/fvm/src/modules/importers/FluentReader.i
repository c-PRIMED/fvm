%module importers
%{
#include "FluentReader.h"
%}

%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"

%include "Mesh.i"
%include "Field.i"

using namespace std;

struct FluentZone 
{
  int ID;
  int iBeg;
  int iEnd;
  int threadType;
  string zoneName;
  string zoneType;
  string zoneVars;
};

struct FluentFaceZone : public FluentZone
{
  int leftCellZoneId;
  int rightCellZoneId;
};

struct FluentCellZone : public FluentZone
{
};

namespace std{
%template(FaceZonesVector) vector<FluentFaceZone*>;
%template(CellZonesVector) vector<FluentCellZone*>;

%template(FaceZonesMap) map<int,FluentFaceZone*>;
%template(CellZonesMap) map<int,FluentCellZone*>;
}

typedef map<int,FluentFaceZone*> FaceZonesMap;
typedef map<int,FluentCellZone*> CellZonesMap;

class FluentReader
{
public:

  FluentReader(const string& fileName);
  void readMesh();
  int getNumCells();
  MeshList getMeshList();
  string getVars();
  
  FaceZonesMap& getFaceZones();
  CellZonesMap& getCellZones();
};
