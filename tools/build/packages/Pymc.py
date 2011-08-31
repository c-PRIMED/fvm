from build_packages import *

class Pymc(BuildPkg):
    requires = ['scipy', 'matplotlib']

    def _installed(self):
        return python_package('pymc', [2,0,0])

    def _build(self):
        import os
        print os.getcwd()
        return self.sys_log("python setup.py config_fc --fcompiler gnu95 build")

    def _install(self):
        do_env("LDFLAGS=")
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

