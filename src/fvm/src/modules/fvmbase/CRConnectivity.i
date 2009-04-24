%{
#include "CRConnectivity.h"

%}

class CRConnectivity
{
public:

 CRConnectivity(const StorageSite& rowSite, const StorageSite& colSite);

  ~CRConnectivity();
 
 int getCount(const int i) const;
	
 int operator()(const int i, const int j) const;

 %extend
{
int __getitem__(int i, int j)
{
return self->operator()(i,j);
}
}

private:
 
 shared_ptr<Array<int> > _row;
 shared_ptr<Array<int> > _col;
};
