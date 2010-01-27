from build_packages import *

class Python(BuildPkg):    
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --enable-shared" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        self.sys_log("make install")
        BuildPkg.pypath = set_python_path(self.blddir)
