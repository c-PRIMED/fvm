from build_packages import *

class Rlog(BuildPkg):
    def _installed(self):
        """
        For some reason not all rlogs are compatible.
        Always use ours for now.

        pathlist=['/usr', '/usr/local', '/opt/local']
        for path in pathlist:
            try:
                f = open(os.path.join(path, 'include', 'rlog', 'rlog.h'), 'r')
            except:
                continue
            if f:
                verbose(2,'Found rlog.h')
                for line in f:
                    if "CURRENT_RLOG_VERSION" in line and '20040503' in line:
                        f.close()
                        return True
                f.close()
        """
        return False

    def _configure(self):
        return self.sys_log("%s/configure --disable-docs --disable-valgrind --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
