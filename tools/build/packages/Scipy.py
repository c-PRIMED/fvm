from build_packages import *

class Scipy(BuildPkg):
    requires = ['numpy', 'python', 'umfpack']

    def _installed(self):
        return python_package('scipy', [0,9])

    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

