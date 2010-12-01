from build_packages import *

# Build GNU time command, because some systems only have BSD time.

class Time(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --program-prefix=g -prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
