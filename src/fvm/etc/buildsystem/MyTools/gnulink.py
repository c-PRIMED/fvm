import SCons.Tool.gnulink
from SCons.Util import CLVar

def generate(env):
    SCons.Tool.gnulink.generate(env)

    env['LINKFLAGS'] = CLVar('-rdynamic')
    if env['DEBUG']:
        env.Append(LINKFLAGS=['-g'])

    env['SHLINKFLAGS'] = env['LINKFLAGS']
    env.Append(SHLINKFLAGS=['-shared'])

    if env['OPENMP']:
        env.Append(SHLINKFLAGS=['-fopenmp'])


    env['SYSLIBS'] = []

def exists(env):
    return SCons.Tool.gnulink.exists(env)

