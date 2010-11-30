from build_packages import *

class Nose(BuildPkg):
    requires = ['python']
    def _installed(self):
        return python_package('nose', [0,10])
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
