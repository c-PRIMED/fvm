from build_packages import *

class Pyodbc(BuildPkg):
    requires = ['python']
    def _configure(self):
        return self.sys_log("python setup.py config")
    def _build(self):
        return self.sys_log("python setup.py build")
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
