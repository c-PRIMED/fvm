from build_packages import *

class Netcdf(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --with-pic -prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        # make -j4 sometimes failed
        return self.sys_log("make")
    def _install(self):
        return self.sys_log("make install")
