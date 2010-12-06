from build_packages import *

class Swig(BuildPkg):
    requires=['numpy']
    def _installed(self):
        try:
            ver = subprocess.Popen('swig -version', \
                                   shell=True, \
                                   stdout=subprocess.PIPE).stdout.read()
            ver  = re.findall(r'SWIG Version ([\d]*)\.*([\d]*)\.*([\d]*)', ver)[0]
            if ver[2] == '': ver[2] = 0
            a,b,c = map(int, ver)
            if a == 2:
                return True
            if a == 1 and b == 3 and c >= 30:
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
