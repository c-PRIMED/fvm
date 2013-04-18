from build_packages import *

class MEMOSA(BuildPkg):
    requires = ['mpi4py', 'h5py', 'nose', 'puq']
    copy_sources = 1
    name = "MEMOSA"
    def _build(self): 
        if sys.platform == 'darwin':
            return False
        hdf5 = self.bld.pkglist['hdf5']
        assert(hdf5)
        hdf5.find_hdf5_vers()
        #return self.sys_log("make -j%s" % (jobs(self.name)))
        return 0
    def _install(self):
        if sys.platform == 'darwin':
            return False
        #os.chdir(self.bdir)
        #self.sys_log("install src/xyz2xdmf %s" % self.bindir)
        os.chdir(self.sdir)
        self.sys_log("install bin/* %s" % self.bindir)
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'memosa'))
        self.sys_log("install -m 644 lib/*.py %s" % os.path.join(self.libdir, 'memosa'))
        self.sys_log("install -m 444 lib/*.h %s" % self.incdir)
        return 0
