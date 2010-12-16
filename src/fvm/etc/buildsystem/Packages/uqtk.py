import distutils.sysconfig


def generate(env):
    env.Append(CPPPATH='$PACKAGESDIR/include')
    env.Append(LIBS=['rlog', 'UQTK', 'netlib-uqtk', 'slatec-uqtk', 
                     'blas-uqtk', 'cvode-uqtk', 'zufall-uqtk', 'lapack'])
    env.Append(LIBPATH=env['PACKAGESLIBDIR'])
