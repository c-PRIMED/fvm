import distutils.sysconfig


def generate(env):
    env.Append(LIBS=['netcdf_c++','netcdf'])

