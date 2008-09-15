#include <string>

using namespace std;

#include "MMReader.h"
#include "AMG.h"

int main()
{
  MMReader reader("/home/sm/tmp/MatrixMarket226.dat",
                  "/home/sm/tmp/rhs226.dat");
  
  shared_ptr<LinearSystem> ls(reader.getLS());

  AMG solver(*ls);
  solver.solve();
  
  return 0;
}
