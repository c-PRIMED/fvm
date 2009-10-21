from build_packages import *

# FLTK (pronounced "fulltick") is a cross-platform C++ GUI toolkit.
# http://www.fltk.org/
# Version 1.1.9 required by gmsh
class Fltk(BuildPkg):
    copy_sources = 1
    def _configure(self):
        return self.sys_log("%s/configure --enable-xft --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
