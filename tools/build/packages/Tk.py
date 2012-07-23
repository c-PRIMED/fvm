from build_packages import *

class Tk(BuildPkg):
    requires=['tcl']    
    def _configure(self):
        if sys.platform == 'darwin':
            dirname = 'macosx'
        else:
            dirname = 'unix'
        tcldir = '%s/../tcl' % (self.sdir)           
        return self.sys_log("%s/%s/configure --with-tcl=%s --prefix=%s --enable-shared" % 
                            (self.sdir, dirname, tcldir, self.blddir))
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        return self.sys_log("make install")
