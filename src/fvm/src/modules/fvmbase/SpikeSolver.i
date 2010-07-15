%{
#include "SpikeSolver.h"
#include "SpikeStorage.h"
%}

class SpikeSolver : public LinearSolver
{
public:

   SpikeSolver(const SpikeStorage& spike_storage);

};

