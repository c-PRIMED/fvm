%{
#include "Model.h"
%}

%include "std_vector.i"


using namespace std;

class Model
{
public:
  Model(const MeshList& meshes);
  virtual void init()=0;
  virtual map<string,shared_ptr<ArrayBase> >&
  getPersistenceData();

  virtual void restart();
};


%template(PersistenceData) std::map<std::string,shared_ptr<ArrayBase> >;
