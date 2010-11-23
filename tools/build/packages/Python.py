from build_packages import *
import config

class Python(BuildPkg):
    def _configure(self):
        if config.config('python', 'Build') == '':
            py = find_executable('python')
            ver = subprocess.Popen('python -V', shell=True, stderr=subprocess.PIPE).stderr.read()
            a, b = re.compile(r'Python ([^.]*).([^.\n]*)').findall(ver)[0]
            if int(a) == 2 and int(b) >= 5:
                self.use_installed = True
                return 'installed'
        self.use_installed = False
        return self.sys_log("%s/configure --prefix=%s --enable-shared" % (self.sdir, self.blddir))
    def _build(self):
        if self.use_installed:
            return "skip"
        return self.sys_log("make")
    def _install(self):
        set_python_path(self.blddir)
        if self.use_installed:
            return "Skip"
        self.sys_log("make install")
