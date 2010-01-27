from build_packages import *
import glob

class Fvm(BuildPkg):
    requires = ['scons', 'rlog', 'swig', 'nose', 'netcdf', 'boost']
    
    def add_required(self, deps):
        x = config('fvm', 'parallel')
        if x != '' and eval(x):
            deps += ['mpi4py', 'parmetis']
        return deps
        
    def status(self):
        x = config('fvm', 'parallel')
        if x != '' and eval(x):
            par = 'parallel'
        else:
            par = 'serial'        
        return '%s-%s' % (par, config('fvm', 'version'))
    
    # from fvm sources
    def get_arch(self):
        if sys.platform == 'linux2':
            if os.uname()[4] == 'ia64':
                return 'lnxia64'
            elif os.uname()[4] == 'x86_64':
                return 'lnx86_64'
            else:
                return 'lnx86'
        elif sys.platform == 'win32':
            return 'ntx86'
        else:
            return sys.platform
    def get_compiler(self, comp):
        vers = '4.2.1'
        if comp == 'intelc':
            comp = 'icc'
            vers = '10.1'
        try:
            line = os.popen("/bin/bash -c '%s --version 2>&1'" % comp).readline()
            vers = comp + '-' + re.compile(r'[^(]*[^)]*\) ([^\n ]*)').findall(line)[0]
        except:
            pass
        return vers
    def _configure(self):            
        pdir = os.path.join(self.sdir, "packages")
        self.sys_log("/bin/mkdir -p %s" % pdir)
        pdir = os.path.join(pdir, self.get_arch())
        return self.sys_log("/bin/ln -fsn %s %s" % (self.blddir, pdir))
    def _build(self):
        par = config(self.name, 'parallel')
        comp = config(self.name, 'compiler')
        ver = config(self.name, 'version')
        self.bld.Scons.set_path(unload=0)
        val = self.sys_log("%s/etc/buildsystem/build -j%s -C %s COMPACTOUTPUT=False PARALLEL=%s VERSION=%s COMPILER=%s" \
                         % (self.sdir, 1 , self.sdir, par, ver, comp))
        self.bld.Scons.set_path(unload=1)
        return val

    def _install(self):
        # The bin directory contains a mixture of libs, python scripts, and binaries.
        # Sort it all out.
        vers = self.get_compiler(config(self.name, 'compiler'))
        rel = config(self.name, 'version')
        pdir = os.path.join(self.sdir, "build", self.get_arch(), vers, rel, "bin")
        os.chdir(pdir)
        all_files = glob.glob('*')
        so_files = glob.glob('*.so')
        py_files = glob.glob('*.py')
        exe_files = [f for f in all_files if f not in so_files and f not in py_files]
        self.sys_log("install %s %s" % (' '.join(exe_files), self.bindir))
        self.sys_log("install *.so %s" % self.libdir)

        # install swig-generated python files in lib/fvm
        self.sys_log("install -d %s" % os.path.join(self.libdir, 'fvm'))
        self.sys_log("install *.py %s " % os.path.join(self.libdir, 'fvm'))
        # install python files from source lib to lib/fvm
        os.chdir(os.path.join(self.sdir, "lib"))
        self.sys_log("install *.py %s " % os.path.join(self.libdir, 'fvm'))

        # install scripts
        os.chdir(os.path.join(self.sdir, "scripts"))
        self.sys_log("install *.py %s" % self.bindir)     
        return 0

    def _clean(self):
        bd = os.path.join(self.sdir, "build")
        pmess("CLEAN", self.name, bd)
        os.system('touch %s/SConstruct' % self.psdir)
        self.pstatus(self.sys_log("rm -rf %s/*" % bd))
