from build_packages import *

def h5py_version():
    cmd = "python -c 'import h5py; print h5py.version.version' 2> /dev/null"
    try:
        ver = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()        
        return re.findall(r'([\d]+)', ver)
    except:
        return []

class H5py(BuildPkg):
    requires = ['numpy', 'hdf5']

    def _installed(self):
        ver = h5py_version()
        if ver and (int(ver[0] > 1)):
            return True
        return False

    def _build(self):
        hdf5 = self.bld.pkglist['hdf5']
        assert(hdf5)
        v, mpi = hdf5.find_hdf5_vers()
        verbose(2, "v=%s mpi=%s" % (v,mpi))
        if v.startswith('1.6'):
            api = '16'
        else:
            api = '18'
        if mpi:
            verbose(1, "Building h5py with mpi")
            do_env('CC=mpicc')
        self.sys_log("python setup.py configure --api=%s" % api)
        ret = self.sys_log("python setup.py build")
        if mpi:
            do_env('CC=mpicc', unload=True)
        return ret

    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
