from build_packages import *

class ParMetis(BuildPkg):
    copy_sources = 1    
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        idir = os.path.join(self.blddir, "include")
        self.sys_log("install --mode=444 parmetis.h %s" % idir)
        return self.sys_log("install *.a *.so *.so.* %s" % self.libdir)
