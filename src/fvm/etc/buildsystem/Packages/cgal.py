import distutils.sysconfig


def generate(env):
    env.Append(LIBS=['CGAL'])

