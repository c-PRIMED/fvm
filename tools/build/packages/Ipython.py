from build_packages import *

# actually a collection of required packages for ipython
class Ipython(BuildPkg):
    requires = ['numpy']
    def _configure(self):
        cfgname = os.path.join(self.bdir,'pkginstall.cfg')
        cf = open(cfgname, 'w')
        cf.write('PREFIX=%s\n' % self.blddir)
        cf.close()
        return 0
    def _build(self):
        return self.sys_log("make")
