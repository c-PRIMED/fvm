from build_packages import *

class Poster(BuildPkg):
    requires = ['python']
    def _install(self):
        #ret = self.sys_log("python setup.py install --prefix=%s" % self.blddir)
        os.chdir(self.sdir)
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'poster'))
        self.sys_log("install -m 644 poster/*.py %s" % os.path.join(self.libdir, 'poster'))
        return 0

