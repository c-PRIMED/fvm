from build_packages import *

class Xdmf(BuildPkg):
    name = 'Xdmf'
    copy_sources = 1
    requires=['numpy']
    def _configure(self):
        return self.sys_log("cmake %s -DCMAKE_INSTALL_PREFIX=%s -DXDMF_WRAP_PYTHON=1 -DCMAKE_BUILD_TYPE=Release" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        os.chdir(os.path.join(self.bdir, "bin"))
        return self.sys_log("install %s *.so %s" % (os.path.join(self.sdir, 'libsrc','Xdmf.py'), self.libdir))
