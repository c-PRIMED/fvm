from build_packages import *

class MEMOSA(BuildPkg):
    requires = ['mpi4py', 'h5py', 'nose', 'puq']
    copy_sources = 1
    name = "MEMOSA"
    def _build(self):
        hdf5 = self.bld.pkglist['hdf5']
        assert(hdf5)
        hdf5.find_hdf5_vers()
        return self.sys_log("make -j%s" % (jobs(self.name)))
    def _install(self):
        os.chdir(self.bdir)
        self.sys_log("install src/xyz2xdmf %s" % self.bindir)
        os.chdir(self.sdir)
        self.sys_log("install bin/* %s" % self.bindir)
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'memosa'))
        self.sys_log("install --mode=644 lib/*.py %s" % os.path.join(self.libdir, 'memosa'))

        idir = os.path.join(self.blddir, "include")
        if not os.path.isdir(idir):
            self.sys_log("/bin/mkdir -p %s" % idir)
        self.sys_log("install --mode=444 lib/*.h %s" % idir)
        return 0
