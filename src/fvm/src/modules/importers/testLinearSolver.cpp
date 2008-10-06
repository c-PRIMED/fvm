#include <string>

using namespace std;

#include "MMReader.h"
#include "AMG.h"

int main(int argc, char *argv[])
{
  MMReader reader(argv[1], argv[2]);
  
  shared_ptr<LinearSystem> ls(reader.getLS());

  AMG solver;
  solver.solve(*ls);
  
  return 0;
}
