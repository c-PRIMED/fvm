from build_packages import *

class Swig(BuildPkg):
    requires=['numpy']
    def _installed(self):
        try:
            ver = subprocess.Popen('swig -version', \
                                   shell=True, \
                                   stdout=subprocess.PIPE).stdout.read()
            a, b = re.findall(r'SWIG Version ([^.]*).([^.\n]*)', ver)[0]
            if int(a) == 2:
                return True
            if int(a) == 1 and int(b) >= 3:
                return True
        except:
            pass
        return False

    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
