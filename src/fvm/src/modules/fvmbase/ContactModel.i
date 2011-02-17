%{
#include "ContactModel.h"
%}


using namespace std;

%import "Model.i"

%include "ContactModel.h"


%template(ContactModelA) ContactModel< ATYPE_STR >;
