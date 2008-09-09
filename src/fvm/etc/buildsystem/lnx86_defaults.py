from SCons.Errors import UserError
from SCons.Util import CLVar

def updateEnv(env):

    #default compiler suite is gcc
    env['CXXVERSION'] = None
    if env['COMPILER'] is None:
        compilerSuite = 'gcc'
    else:
        c = env['COMPILER'].split('-')
        compilerSuite=c[0]
        if len(c) > 1:
            env['CXXVERSION'] = c[1]
        
    if compilerSuite == 'gcc':
        env.loadTools( [ 'gnulink', 'g++'] )
    elif compilerSuite == 'intelc':
        env.loadTools(['gnulink', 'intelc'])
    else:
        raise UserError, 'Unknown compiler: %s' % env['COMPILER']

    env.AppendUnique(CPPDEFINES = CLVar('OS_LINUX'))
