%{
#include "LinearSystemMerger.h"
%}

using namespace std;

class LinearSystemMerger
{
public:

   LinearSystemMerger( int target_proc_id, const set<int>& group, LinearSystem& ls );
   
   void debug_print();


};
