from build_packages import *

class Umfpack(BuildPkg):
    requires = ['amd']
    def _build(self):
        return self.sys_log("INSTALL_LIB=%s INSTALL_INCLUDE=%s make -j%s library" % (self.libdir, self.incdir, jobs(self.name)))
    def _install(self):
        return self.sys_log("INSTALL_LIB=%s INSTALL_INCLUDE=%s make install" % (self.libdir, self.incdir))
