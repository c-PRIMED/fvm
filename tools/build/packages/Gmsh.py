from build_packages import *

class Gmsh(BuildPkg):
    requires = ['fltk', 'fftw']
    copy_sources = 1    
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --with-fltk-prefix= %s --with-gsl-prefix=%s" \
                                % (self.sdir, self.blddir, self.blddir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
