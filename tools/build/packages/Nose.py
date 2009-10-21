from build_packages import *

class Nose(BuildPkg):
    requires = ['python']
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
