from build_packages import *

class Puq(BuildPkg):
    requires = ['mpi4py', 'h5py', 'scipy', 'matplotlib', 'nose', 'jsonpickle', 'sympy', 'time']
    copy_sources = 1
    def _build(self):
        return self.sys_log("make -j%s" % (jobs(self.name)))
    def _install(self):
        os.chdir(os.path.join(self.bdir, 'adap', 'src'))
        self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        os.chdir(self.bdir)
        self.sys_log("install src/sparse_grid_cc %s" % self.bindir)
        os.chdir(self.sdir)
        self.sys_log("install bin/* %s" % self.bindir)
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'puq'))
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'puqutil'))
        self.sys_log("install -m 644 lib/*.py %s" % os.path.join(self.libdir, 'puq'))
        self.sys_log("install -m 644 puqutil/*.py %s" % os.path.join(self.libdir, 'puqutil'))

        idir = os.path.join(self.blddir, "include")
        if not os.path.isdir(idir):
            self.sys_log("/bin/mkdir -p %s" % idir)
        self.sys_log("install -m 444 puqutil/*.h %s" % idir)
        return 0
