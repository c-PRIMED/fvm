import os, re
import SCons.Tool
from SCons.Util import CLVar
from Arch import is64Bit

__import__('SCons.Tool', globals(), locals(), ['c++'])
cppTool = getattr(SCons.Tool, 'c++')

def generate(env):
    cppTool.generate(env)

    if not env.get('CXXVERSION'):
        try:
            line = os.popen("/bin/bash -c 'gcc --version 2>&1'").readline()
            env['CXXVERSION'] = line.split()[2]
        except:
            env['CXXVERSION'] = '4.2.1'

    env['COMPILER'] = 'gcc-' + env['CXXVERSION']

    env['CXXFLAGS'] = CLVar('-Wall -Woverloaded-virtual -ftemplate-depth-200')

    if is64Bit():
        env.Append(CPPDEFINES=CLVar('OS_64BIT'))

    if env['DEBUG']:
        env.Append(CXXFLAGS=['-g', '-O0'])
    else:
        env.Append(CXXFLAGS=['-O3', '-finline-limit=500'])


    if env['OPENMP']:
        env.Append(CXXFLAGS=['-fopenmp'])
    
    env['CCFLAGS'] = env['CXXFLAGS']

    env['SHCXXFLAGS'] = CLVar('$CXXFLAGS -fPIC')

def exists(env):
    return env.Detect(['g++'])

