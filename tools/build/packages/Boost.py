from build_packages import *

# We don't build boost, just use the headers
class Boost(BuildPkg):
    def _install(self):
        idir = os.path.join(self.blddir, "include")
        return self.sys_log("/bin/ln -fs %s %s" % (self.bdir, idir))


