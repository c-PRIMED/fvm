from build_packages import *

class Openpyxl(BuildPkg):
    requires = ['python', 'setuptools']
    def _installed(self):
        return python_package('openpyxl', [1,6])
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
