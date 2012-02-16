from build_packages import *

# Build GMP, The GNU Multiple Precision Library http://gmplib.org/

class Gmp(BuildPkg):
    def _installed(self):
        pathlist=['/usr', '/usr/local', '/opt/local']
        for path in pathlist:
            try:
                f = open(os.path.join(path, 'include', 'gmp.h'), 'r')
            except:
                continue
            if f:
                f.close()                
                verbose(2,'Found gmp.h')
                return True
        return False

    def _configure(self):
        return self.sys_log("%s/configure -prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
