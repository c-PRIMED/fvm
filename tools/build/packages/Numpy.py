from build_packages import *

class Numpy(BuildPkg):
    requires = ['python']

    def _installed(self):
        return python_package('numpy', [1,6,2])
    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

