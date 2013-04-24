from build_packages import *

class Puq(BuildPkg):
    requires = ['setuptools', 'h5py', 'scipy', 'matplotlib', 'nose', 'jsonpickle', 'sympy', 'time', 'swig', 'poster', 'pymc']

    def _install(self):
        os.chdir(self.sdir)        
        idir = os.path.join(self.blddir, "include")
        if not os.path.isdir(idir):
            self.sys_log("/bin/mkdir -p %s" % idir)
        self.sys_log("install -m 444 puqutil/*.h %s" % idir)
        self.sys_log("python setup.py clean") 
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)
