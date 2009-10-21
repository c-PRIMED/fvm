from build_packages import *

class MEMOSA(BuildPkg):
    requires=['mpi4py']
    copy_sources = 1    
    name = "MEMOSA"
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        os.chdir(self.bdir)
        self.sys_log("install src/sparse_grid_cc %s" % self.bindir)
        os.chdir(self.sdir)
        self.sys_log("install bin/* %s" % self.bindir)
        self.sys_log("install lib/* %s" % self.libdir)        
        return 0
