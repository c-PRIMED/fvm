from build_packages import *

class Pytables(BuildPkg):
    requires = ['hdf5', 'numpy']
    def _install(self):
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        return ret
