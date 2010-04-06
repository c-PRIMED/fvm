from build_packages import *

class Vitables(BuildPkg):
    requires = ['pytables']
    def _install(self):
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        return ret

