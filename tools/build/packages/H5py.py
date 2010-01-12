from build_packages import *


class H5py(BuildPkg):
    requires = ['numpy', 'hdf5']
    def find_h5_inc(self):
        f = ''
        pathlist = ['/usr/include', '/usr/local/include']
        if os.environ.has_key('HDF5_DIR'):
            pathlist.append(os.path.join(os.environ['HDF5_DIR'], 'include'))
        for path in pathlist:
            verbose(2,'Checking for HDF5 headers in %s' % path)
            for fn in ['H5pubconf.h', 'H5pubconf-64.h', 'H5pubconf-32.h']:
                try:
                    f = open(os.path.join(path, fn), 'r')
                except:
                    pass
                if f:
                    verbose(2,'Found HDF5 headers.')
                    return f
        if not f:
            fatal("\nhdf5 include files not found. Please install them and restart build.")
        return f

    def h5_vers(self, f):
        if not f:
            return ''
        v = ''
        mpi = ''
        while True:
            line = f.readline()
            if not line:
                break
            if not v:
                v = re.findall(r'#define H5_PACKAGE_VERSION "([^"]*)"', line)
                if v:
                    v = v[0]
            if not mpi:
                mpi = re.findall(r'#define H5_HAVE_MPI_GET_SIZE ([0-9]*)', line)
            if v and mpi:
                break
        verbose(2, 'HDF5 version=%s, mpi=%s' % (v, mpi))
        return v, mpi

    def _build(self):
        v = config(self.name, 'HDF5_VERSION')
        if not v:
            v, mpi = self.h5_vers(self.find_h5_inc())
        if not v:
            warning("HDF5_VERSION not set and hdf5 include files not found.")
            warning("Assuming hdf5 version = 1.8 and continuing...")
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
