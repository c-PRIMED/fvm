import distutils.sysconfig


def generate(env):
    env.AppendUnique(CPPPATH='$PACKAGESDIR/include')
    env.Append(LIBS=['parmetis','metis','mpi_cxx','mpi'])
    env.AppendUnique(LIBPATH=env['PACKAGESLIBDIR'])
