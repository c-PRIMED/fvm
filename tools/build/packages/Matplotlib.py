from build_packages import *

class Matplotlib(BuildPkg):
    requires = ['numpy', 'freetype']

    def _installed(self):
        # need 1.0.1
        if python_package('matplotlib', [1,0,1]):
            return True
        return False

    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

