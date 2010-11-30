from build_packages import *

class Pytables(BuildPkg):
    requires = ['h5py', 'numpy']
    def _install(self):
        hdf5 = self.bld.pkglist['hdf5']
        assert(hdf5)
        hdf5.find_hdf5_vers()
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        return ret
