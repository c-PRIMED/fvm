from build_packages import *

class Numpy(BuildPkg):
    requires = ['python']

    def _installed(self):
        # We need our modified numpy on MacOSX
        if sys.platform == 'darwin':
            return False
        return python_package('numpy', [1,3])
    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

