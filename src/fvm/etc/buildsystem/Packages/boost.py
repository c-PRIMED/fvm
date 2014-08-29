import distutils.sysconfig
import os

def generate(env):
    print 'Package BOOST'
    try:
        cgi = os.environ['BOOST_INCLUDE']
        if cgi.startswith('-I'):
            cgi = cgi.lstrip('-I')

        if cgi:
            env.Append(CPPPATH=[cgi])
    except:
        pass

    #env.Append(LIBS=['BOOST'])


