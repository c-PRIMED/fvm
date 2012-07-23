from build_packages import *

class Python(BuildPkg):
    requires=['tk']    
    def _installed(self):
        ver = subprocess.Popen('python -V', shell=True, stderr=subprocess.PIPE).stderr.read()
        a, b = re.findall(r'Python ([^.]*).([^.\n]*)', ver)[0]
        if int(a) == 2 and int(b) >= 6:
            return True
        return False
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --enable-shared" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        ret = self.sys_log("make install")
        set_python_path(self.blddir, reset=True)
        return ret
