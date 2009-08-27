import distutils.sysconfig


def generate(env):
    env.AppendUnique(CPPPATH='$PACKAGESDIR/include')
    env.Append(LIBS=['parmetis','metis'])
    env.AppendUnique(LIBPATH=env['PACKAGESLIBDIR'])
