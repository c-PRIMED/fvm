from build_packages import *

class Tcl(BuildPkg):
    def _configure(self):
        if sys.platform == 'darwin':
            dirname = 'macosx'
        else:
            dirname = 'unix'        
        return self.sys_log("%s/%s/configure --prefix=%s --enable-shared" % 
                            (self.sdir, dirname, self.blddir))
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        return self.sys_log("make install")
