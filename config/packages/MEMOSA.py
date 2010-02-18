from build_packages import *

class MEMOSA(BuildPkg):
    requires = ['mpi4py', 'h5py', 'scipy', 'matplotlib', 'nose']
    copy_sources = 1    
    name = "MEMOSA"
    def _build(self):
        return self.sys_log("HDF5_DIR=%s make -j%s" % (self.blddir, jobs(self.name)))
    def _install(self):
        os.chdir(self.bdir)
        self.sys_log("install src/sparse_grid_cc %s" % self.bindir)
        self.sys_log("install src/xyz2xdmf %s" % self.bindir)
        os.chdir(self.sdir)
        self.sys_log("install bin/* %s" % self.bindir)
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'memosa'))
        self.sys_log("install --mode=644 lib/* %s" % os.path.join(self.libdir, 'memosa'))
        return 0
