from build_packages import *

class Netcdf(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --with-pic -prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
