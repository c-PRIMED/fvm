%{
#include "StorageSiteMerger.h"
%}

using namespace std;

class StorageSiteMerger
{
public:

   StorageSiteMerger( int target_proc_id, const set<int>& group, const StorageSite& site );
   
   void debug_print();

   int  getSelfCount()  const ;
   int  getGhostCount() const ;
   int  getCount()      const ;

};
