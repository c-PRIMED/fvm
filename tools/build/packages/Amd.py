from build_packages import *

class Amd(BuildPkg):
    requires = ['ufconfig']
    def _install(self):
        # Need to rename the unpack directory because UMFPACK
        # is hardcoded to look for it with a certain name
        newdir = os.path.join(os.path.split(self.sdir)[0], 'AMD')
        shutil.move(self.sdir, newdir)
        self.sdir = newdir
        self.bdir = newdir
        return self.sys_log("INSTALL_LIB=%s INSTALL_INCLUDE=%s make install" % (self.libdir, self.incdir))
    def _build(self):
        return self.sys_log("INSTALL_LIB=%s INSTALL_INCLUDE=%s make" % (self.libdir, self.incdir))
