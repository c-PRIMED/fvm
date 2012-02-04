from build_packages import *
import config

class Flann(BuildPkg):
    requires = ['hdf5']

    def _installed(self):
        return find_lib('flann', '1.7.1')

    def _configure(self):
        return self.sys_log("cmake %s -DBUILD_CUDA_LIB=OFF -DBUILD_MATLAB_BINDINGS=OFF -DCMAKE_INSTALL_PREFIX=%s -DCMAKE_BUILD_TYPE=Release -DLATEX_OUTPUT_PATH=/tmp" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")


