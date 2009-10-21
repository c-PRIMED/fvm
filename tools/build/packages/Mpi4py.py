from build_packages import *

class Mpi4py(BuildPkg):
    requires = ['numpy']

    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)

