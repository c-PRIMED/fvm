from build_packages import *

class Matplotlib(BuildPkg):
    requires = ['numpy']

    def _installed(self):
        # 0.99 is good. 1.0 has a problem when plotting
        if python_package('matplotlib', [0,99]) \
            and not python_package('matplotlib', [1,0]):
            return True
        return False

    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

