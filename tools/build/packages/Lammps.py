from build_packages import *

class Lammps(BuildPkg):
    requires = ['fftw']
    copy_sources = 1    
    def _build(self):
        os.chdir('src')
        ret = self.sys_log("make -j%s" % jobs(self.name));
        os.chdir('..')
        return ret
    def _install(self):
        os.chdir('src')
        self.sys_log("install lmp_* %s" % self.bindir)
        # for testing, we want a consistent name for the executable
        ret = self.sys_log("/bin/ln -fs %s/lmp_* %s/lammps" % (self.bindir, self.bindir))
        os.chdir('..')
        return ret
