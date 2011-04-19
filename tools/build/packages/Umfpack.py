from build_packages import *

class Umfpack(BuildPkg):
    requires = ['amd']
    def _build(self):
        idir = os.path.join(self.blddir, "include")
        return self.sys_log("INSTALL_LIB=%s INSTALL_INCLUDE=%s make -j%s library" % (self.blddir, idir, jobs(self.name)))
    def _install(self):
        idir = os.path.join(self.blddir, "include")
        return self.sys_log("INSTALL_LIB=%s INSTALL_INCLUDE=%s make install" % (self.blddir, idir))
