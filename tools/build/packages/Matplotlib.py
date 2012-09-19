from build_packages import *

class Matplotlib(BuildPkg):
    requires = ['numpy', 'freetype']

    def _installed(self):
        # need 1.0.1
        if python_package('matplotlib', [1,1,0]):
            return True
        return False

    def _install(self):
        do_env("LDFLAGS=")
        if sys.platform != 'darwin':
            cfgname = os.path.join(self.sdir, 'setup.cfg')
            f = open(cfgname, 'w')
            print >>f, "[directories]"
            print >>f, "basedirlist = /usr /usr/local %s" % self.blddir
            f.close()
        ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        do_env("LDFLAGS=", unload=True)
        return ret

