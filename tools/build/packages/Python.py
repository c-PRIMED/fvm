from build_packages import *

class Python(BuildPkg):
    def _installed(self):
        ver = subprocess.Popen('python -V', shell=True, stderr=subprocess.PIPE).stderr.read()
        a, b = re.findall(r'Python ([^.]*).([^.\n]*)', ver)[0]
        if int(a) == 2 and int(b) >= 5:
            return True
        return False
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --enable-shared" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        set_python_path(self.blddir)
        self.sys_log("make install")
