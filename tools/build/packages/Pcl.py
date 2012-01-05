from build_packages import *
import config

class Pcl(BuildPkg):
    requires = ['flann', 'eigen', 'boost']
    
    def _configure(self):
        return self.sys_log("cmake %s -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_BUILD_TYPE=Release" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")


