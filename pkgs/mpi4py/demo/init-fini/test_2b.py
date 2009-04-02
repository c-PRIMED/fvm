from mpi4py import rc
rc.initialize = False

from mpi4py import MPI
assert not MPI.Is_initialized()
assert not MPI.Is_finalized()

MPI.Init_thread()
assert MPI.Is_initialized()
assert not MPI.Is_finalized()

name, _ = MPI.get_vendor()
if name == 'MPICH2':
    assert MPI.Query_thread() == MPI.THREAD_MULTIPLE

MPI.Finalize()
assert MPI.Is_initialized()
assert MPI.Is_finalized()
