from build_packages import *

class ParMetis(BuildPkg):
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        self.sys_log("install -m 444 parmetis.h %s" % self.incdir)
        return self.sys_log("install *.a *.so* %s" % self.libdir)
