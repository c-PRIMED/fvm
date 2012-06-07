from build_packages import *

class Jsonpickle(BuildPkg):
    requires = ['python']
    def _install(self):
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        return ret

