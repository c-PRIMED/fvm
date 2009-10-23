from build_packages import *

# We don't build boost, just use the headers
class Boost(BuildPkg):
    pass
# What was done below is crazy. If boost is needed, the proper
# solution is to implement copy_srcs() method for this package
# to directly untar the headers into the include directory.
# See MPM.
#    def _install(self):
#        idir = os.path.join(self.blddir, "include")
#        return self.sys_log("/bin/ln -fs %s %s" % (self.bdir, idir))


