// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#include <string>

using namespace std;

#include "MMReader.h"
#include "AMG.h"
#ifdef FVM_PARALLEL
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{
#ifdef FVM_PARALLEL
  MPI::Init(argc, argv);
#endif
  MMReader reader(argv[1], argv[2]);
  
  shared_ptr<LinearSystem> ls(reader.getLS());

  AMG solver;
  solver.solve(*ls);
  
  return 0;
}
