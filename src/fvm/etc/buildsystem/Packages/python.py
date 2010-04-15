from distutils.sysconfig import get_config_vars, get_python_version, get_python_inc
import numpy, os.path
import Arch
def generate(env):
    cv = get_config_vars()
    pylib = 'python'+get_python_version()
    pylibpath = [cv[var] for var in ('LIBDIR', 'LIBPL')]
    env.Append(CPPPATH=[get_python_inc(),numpy.get_include()])
    if Arch.isWindows():
        env.Append(LIBS=pylib.replace(".",""))
        env.Append(LIBPATH=pylibpath)
    else:
        env.Append(LIBS=pylib)
        env.Append(LIBPATH=pylibpath)
