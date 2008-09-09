import distutils.sysconfig


def generate(env):
    env.Append(CPPPATH='$PACKAGESDIR/include')
    env.Append(LIBS=['rlog'])
    env.Append(LIBPATH=env['PACKAGESLIBDIR'])
