from build_packages import *

class Fvm(BuildPkg):
    requires = ['scons', 'rlog', 'boost', 'swig', 'nose', 'netcdf']
    
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
    def getArch(self):
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
    def getCompiler(self, comp):
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
        pdir = os.path.join(pdir, self.getArch())
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
        vers = self.getCompiler(config(self.name, 'compiler'))
        rel = config(self.name, 'version')
        pdir = os.path.join(self.sdir, "build", self.getArch(), vers, rel, "bin")
        os.chdir(pdir)
        self.sys_log("install testLinearSolver %s" % self.bindir)
        self.sys_log("install *.py *.so %s" % self.libdir)        
        # install scripts
        pdir = os.path.join(self.sdir, "scripts")
        os.chdir(pdir)   
        self.sys_log("install *.py %s" % self.bindir)     
        return 0

    def _clean(self):
        bd = os.path.join(self.sdir, "build")
        pmess("CLEAN", self.name, bd)
        os.system('touch %s/SConstruct' % self.psdir)
        return self.sys_log("rm -rf %s/*" % bd)
