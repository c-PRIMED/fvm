import distutils.sysconfig
import os

def generate(env):
    try:
        cgi = os.environ['CGAL_INCLUDE']
        if cgi.startswith('-I'):
            cgi = cgi.lstrip('-I')

        cgd = os.environ['CGAL_DIR']
        if cgd.endswith('/CGAL'):
            cgd = cgd.rstrip('/CGAL')

        if cgd:
            env.Append(LIBPATH=cgd)
        
        if cgi:
            env.Append(CPPPATH=[cgi])
    except:
        pass

    env.Append(LIBS=['CGAL'])
    if env['DEBUG']:
        env.Append(CXXFLAGS=['-DCGAL_DISABLE_ROUNDING_MATH_CHECK'])


