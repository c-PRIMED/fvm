from build_packages import *
import config

class Hdf5(BuildPkg):

    def find_hdf5_inc(self, installed=False, in_build=False):
        f = ''
        if installed:
            pathlist=['/usr', '/usr/local']
            if os.environ.has_key('HDF5_DIR'):
                pathlist.insert(0, os.environ['HDF5_DIR'])
        elif in_build:
            pathlist=[self.blddir]
        else:
            pathlist=[self.blddir, '/usr', '/usr/local']
            if os.environ.has_key('HDF5_DIR'):
                pathlist.insert(0, os.environ['HDF5_DIR'])

        for path in pathlist:
            verbose(2,'Checking for HDF5 headers in %s' % path)
            for fn in ['H5pubconf.h', 'H5pubconf-64.h', 'H5pubconf-32.h']:
                try:
                    f = open(os.path.join(path, 'include', fn), 'r')
                except:
                    pass
                if f:
                    verbose(2,'Found HDF5 headers.')
                    if path != '/usr':
                        do_env('HDF5_DIR=%s' % path)
                    return f
        return f

    def find_hdf5_vers(self, installed=False, in_build=False):
        v = ''
        mpi = ''
        f = self.find_hdf5_inc(installed, in_build)
        if f:
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

    def _configure(self):
        # Do we really need to build HDF5 or can we use the installed one?
        print os.environ['CC']
        print os.environ['CXX']
        print os.environ['FC']
        
        x = config.config(self.name, 'Build')
        if x != '' and eval(x):
            verbose(2, "HDF5 MUST BUILD")
        else:
            v, mpi = self.find_hdf5_vers(in_build=True)
            if not v:
                v, mpi = self.find_hdf5_vers(installed=True)
                if v:
                    self.use_installed = True
                    return "installed"
        self.use_installed = False
        return self.sys_log("%s/configure --enable-fortran --enable-cxx --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        if self.use_installed:
            return "skip"
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        if self.use_installed:
            return "skip"
        return self.sys_log("make install")


