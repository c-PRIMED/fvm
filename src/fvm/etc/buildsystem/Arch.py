import os
import sys

def getArch():
    if sys.platform == 'linux2':
        if os.uname()[4] == 'ia64':
            return 'lnxia64'
        elif os.uname()[4] == 'x86_64':
            return 'lnx86_64'
        else:
            return 'lnx86'
    elif sys.platform == 'win32':
        return 'ntx86'
    else:
        return sys.platform

def is64Bit():
    return getArch() in ['lnx86_64']

def isWindows():
    return sys.platform == 'win32'

def isLinux():
    return sys.platform in ['lnx86', 'lnx86_64']

def updateEnv(env):

    if not env.get('ARCH'):
        env['ARCH'] = getArch()

    try:
        archDefaults = __import__(env['ARCH']+'_defaults')
    except ImportError:
        raise ImportError('Defaults for "%s not found"' % env['ARCH'])

    archDefaults.updateEnv(env)

    if env['COMPILER'] is None:
        env['COMPILER'] = "%s-%s" % (env['CXXTOOL'], env['CXXVERSION'])

