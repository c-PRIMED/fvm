from build_packages import *

class Rlog(BuildPkg):
    copy_sources = 1    
    def _configure(self):
        return self.sys_log("%s/configure --disable-docs --disable-valgrind --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
