from build_packages import *

class Setuptools(BuildPkg):
    requires = ['python']
    def _installed(self):
        return python_package('setuptools', [0,6])
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
