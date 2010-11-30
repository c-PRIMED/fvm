from build_packages import *

class Matplotlib(BuildPkg):
    requires = ['numpy']

    def _installed(self):
        return python_package('matplotlib', [0,99])

    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

