from distutils.sysconfig import get_config_vars, get_python_version, get_python_inc
import os.path
import Arch
def generate(env):
    cv = get_config_vars()
    pylib = 'python'+get_python_version()
    env.Append(CPPPATH=get_python_inc())
    if Arch.isWindows():
        env.Append(LIBS=[pylib.replace(".","")])
        env.Append(LIBPATH=os.path.join(cv['prefix'],"libs"))
    else:
        env.Append(LIBS=[pylib])
        env.Append(LIBPATH=env['PACKAGESLIBDIR'])
