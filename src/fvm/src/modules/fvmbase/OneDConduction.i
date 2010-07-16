#ifndef _ONEDCONDUCTION_H_
#define _ONEDCONDUCTION_H_

%{
#include "OneDConduction.h"
%}

using namespace std;

template<class T>
class OneDConduction
{
public:
  OneDConduction(const int nCells, const T& kConst)
  {}
  
  void solve();
    %extend
  {
    boost::shared_ptr<ArrayBase> getSolution()
    {
      return self->getSolution();
    }
  }
};


%template(OneDConductionA) OneDConduction< ATYPE_STR >;

#endif
