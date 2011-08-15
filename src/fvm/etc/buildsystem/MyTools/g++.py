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
            env['CXXVERSION'] = re.compile(r'[^(]*[^)]*\) ([^\n ]*)').findall(line)[0]
        except:
            env['CXXVERSION'] = '4.2.1'

    env['COMPILER'] = 'gcc-' + env['CXXVERSION']

    env['CXXFLAGS'] = CLVar('-Wall -fno-strict-aliasing -Woverloaded-virtual -ftemplate-depth-200 -frounding-math')

    if is64Bit():
        env.Append(CPPDEFINES=CLVar('OS_64BIT'))

    if env['DEBUG']:
        env.Append(CXXFLAGS=['-g', '-O0'])
    elif env['PROFILE']:
        env.Append(CXXFLAGS=['-fprofile-arcs', '-ftest-coverage'])
    else:
        env.Append(CXXFLAGS=['-O3', '-finline-limit=500'])


    if env['OPENMP']:
        env.Append(CXXFLAGS=['-fopenmp'])

    if env['PARALLEL']:
        env['CXX'] = 'mpicxx'
        env.Append(CXXFLAGS=['-DFVM_PARALLEL'])
        #bug fix for mpich
        env.Append(CXXFLAGS=['-DMPICH_IGNORE_CXX_SEEK'])
    else:
        env['CXX'] = 'c++'

    #env['CCFLAGS'] = env['CXXFLAGS']

    env['SHCXXFLAGS'] = CLVar('$CXXFLAGS -fPIC')

def exists(env):
    if env['PARALLEL']:
        return env.Detect(['mpicxx'])
    else:
        return env.Detect(['g++'])

