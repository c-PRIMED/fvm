# Build packages
# Martin Hunt <mmh@purdue.edu>

"""
Build package definitions.
"""

import sys, os, testing, cgi, string, re
from build_utils import *
from config import config

def create_build_dir():
    for p in [BuildPkg.blddir,
              BuildPkg.logdir,
              BuildPkg.libdir,
              BuildPkg.bindir,
              ]:

        if not os.path.isdir(p):
            try:
                os.makedirs(p)
            except:
                fatal("error creating directory " + p)

# abstract superclass
class BuildPkg:

    def setup(cname, srcpath):
        BuildPkg.topdir = srcpath
        BuildPkg.blddir = os.path.join(os.getcwd(), "build-%s" % cname)
        BuildPkg.logdir = os.path.join(BuildPkg.blddir, "log")
        BuildPkg.bindir = os.path.join(BuildPkg.blddir, "bin")
        BuildPkg.libdir = os.path.join(BuildPkg.blddir, "lib")
        create_build_dir()
        
        BuildPkg.packages = [
            Python("pkgs/python", 0),
            Numpy("pkgs/numpy", 1),
            Ipython("pkgs/ipython.tgz",-1),
            Mpi4py("pkgs/mpi4py-1.1.0.bz2", -1),
            Nose("pkgs/python-nose", 1),            
            Gsl("pkgs/gsl", 0),            
            Fltk("pkgs/fltk", 1),
            Gmsh("pkgs/gmsh", 1),
            Rlog("pkgs/rlog", 1),
            Boost("pkgs/boost.tgz", -1),
            Swig("pkgs/swig", 0),
            Fftw("pkgs/fftw", 0),
            H5py("pkgs/h5py-1.2.0.bz2", -1),
            Xdmf("pkgs/Xdmf", 0),
            Netcdf("pkgs/netcdf", 2),
            NetCDF4("pkgs/netCDF4-0.8.1.bz2", -1),
            ParMetis("pkgs/ParMetis", 1),
            Lammps("src/lammps", 1),
            MPM("src/MPM", 2),
            Fvm("src/fvm", 0),
            MEMOSA("src/MEMOSA", 0),
            ]
        
    setup = staticmethod(setup)

    def __init__(self, sdir, copytype):
        if not hasattr(self, 'name'):
            self.name = string.lower(self.__class__.__name__)
        self.sdir = os.path.join(self.topdir, sdir)
        self.copy_sources = copytype
        self.bdir = os.path.join(self.blddir, "build", self.name)

    def pstatus(self, state, option=''):
        "update build status"
        if state == 0 or state == None:
            cprint('GREEN', 'ok ' + option)
        elif isinstance(state, int):
            cprint('RED', 'failed ' + option)
            print "Check contents of", self.logfile
            print "\n"
            os.system("tail -40 " + self.logfile)
            sys.exit(state)
        else:
            cprint('YELLOW', state)

    def configure(self):
        self.state = 'configure'
        self.logfile = os.path.join(self.logdir, self.name+"-conf.log")
        remove_file(self.logfile)
        pmess("CONF",self.name,self.bdir)
        # remove any old sources
        os.system("/bin/rm -rf %s" % self.bdir)            
        # get new sources
        copytree(self.sdir, self.bdir, self.copy_sources)
        self.sdir = self.sdir.split('.')[0]
        os.chdir(self.bdir)
        run_commands(self.name, 'before')
        self.pstatus(self._configure())

    def clean(self):
        self.state = 'clean'
        self.logfile = ''
        os.chdir(self.bdir)
        self._clean()
        if self.sdir != self.bdir:
            os.chdir(self.sdir)
            if os.path.isdir(self.bdir):            
                os.system("/bin/rm -rf %s/*" % self.bdir)

    def build(self):
        self.state = 'build'
        self.logfile = os.path.join(self.logdir, self.name+"-build.log")
        remove_file(self.logfile)
        pmess("BUILD",self.name,self.bdir)
        os.chdir(self.bdir)
        self.pstatus(self._build())

    def install(self):
        self.state = 'install'
        self.logfile = os.path.join(self.logdir, self.name+"-install.log")
        remove_file(self.logfile)
        pmess("INSTALL",self.name,self.blddir)
        os.chdir(self.bdir)
        self.pstatus(self._install())
        run_commands(self.name, 'after')

    def test(self):
        self.state = 'testing'
        self.logfile = os.path.join(self.logdir, self.name+"-test.xml")
        remove_file(self.logfile)
        pmess("TEST",self.name,self.blddir)
        ok, errs = self._test()
        if errs:
            cprint('YELLOW', "%s OK, %s FAIL" % (ok, errs))
        else:
            cprint('GREEN', "%s OK" % ok)
        run_commands(self.name,'after')
        return ok, errs
        
    def sys_log(self, cmd):
        "Execute a system call and log the result."
        # get configuration variable
        e = config(self.name,self.state)
        e = e.replace('BUILDDIR', self.blddir)
        e = e.replace('SRCDIR', self.sdir)
        e = e.replace('TMPBDIR', self.bdir)
        e = e.replace('LOGDIR', self.logdir)
        cmd = cmd + " " + e
        debug(cmd)
        if self.logfile != '':
            f = open(self.logfile, 'a')
            print >> f,"EXECUTING:", cmd
            f.close()
            cmd = cmd + " >>" + self.logfile + " 2>&1"
        return os.system("/bin/bash -c '%s'" % cmd)

    # subclasses must redefine these
    def _configure(self):
        pass
    def _clean(self):
        pass
    def _build(self):
        pass
    def _install(self):
        pass
    def _test(self):
        return testing.do_tests(self.name, self.sdir, self.logfile)
        

#########################################################
# PACKAGES 
#
# These contain the build methods for each package.
# These will mostly be the same, but we implement
# them this way to allow maximum flexibility.
#########################################################

# actually a collection of required packages for ipython
class Ipython(BuildPkg):
    def _configure(self):
        cfgname = os.path.join(self.bdir,'pkginstall.cfg')
        cf = open(cfgname, 'w')
        cf.write('PREFIX=%s\n' % self.blddir)
        cf.close()
        return 0
    def _build(self):
        return self.sys_log("make");
        
class Lammps(BuildPkg):
    def _build(self):
        os.chdir('src')
        ret = self.sys_log("make -j%s" % jobs(self.name));
        os.chdir('..')
        return ret
    def _clean(self):
        os.chdir('src')
        self.sys_log("make clean-all")
        self.sys_log("/bin/rm -f lmp_*")
        os.chdir('..')
        return 0
    def _install(self):
        os.chdir('src')
        self.sys_log("install lmp_* %s" % self.bindir)
        # for testing, we want a consistent name for the executable
        ret = self.sys_log("/bin/ln -fs %s/lmp_* %s/lammps" % (self.bindir, self.bindir))
        os.chdir('..')
        return ret

class Gmsh(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --with-fltk-prefix= %s --with-gsl-prefix=%s" \
                                % (self.sdir, self.blddir, self.blddir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class H5py(BuildPkg):
    def find_h5_inc(self):
        f = ''
        for path in ['/usr/include','/usr/local/include']:
            for fn in ['H5pubconf.h', '/usr/include/H5pubconf-64.h']:
                try:
                    f = open(os.path.join(path, fn), 'r')
                except:
                    pass
                if f:
                    return f
        if not f:
            fatal("\nhdf5 include files not found. Please install them and restart build.")
        return f

    def h5_vers(self,f):
        if not f:
            return ''
        v = ''
        mpi=''
        while True:
            line = f.readline()
            if not line:
                break
            if not v:
                v = re.findall(r'#define H5_PACKAGE_VERSION "([^"]*)"',line)
                if v:
                    v = v[0]
            if not mpi:
                mpi = re.findall(r'#define H5_HAVE_MPI_GET_SIZE ([0-9]*)',line)
            if v and mpi:
                break
        return v,mpi
    
    def _build(self):
        v = config(self.name,'HDF5_VERSION')
        if not v:
            v, mpi = self.h5_vers(self.find_h5_inc())
        if not v:
            warning("HDF5_VERSION not set and hdf5 include files not found.")
            warning("Assuming hdf5 version = 1.8 and continuing...")
        if v.startswith('1.6'):
            api = '16'
        else:
            api = '18'
        if mpi:
            verbose("Building h5py with mpi")
            do_env('CC=mpicc')
        self.sys_log("python setup.py configure --api=%s" % api)
        ret = self.sys_log("python setup.py build")
        if mpi:
            do_env('CC=mpicc', unload=True)
        return ret
    
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)

class Numpy(BuildPkg):
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)

class Nose(BuildPkg):
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)

class Mpi4py(BuildPkg):
    def _install(self):
        return self.sys_log("python setup.py install --prefix=%s" % self.blddir)

class ParMetis(BuildPkg):
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        idir = os.path.join(self.blddir, "include")
        self.sys_log("install --mode=444 parmetis.h %s" % idir)
        return self.sys_log("install *.a *.so *.so.* %s" % self.libdir)

class NetCDF4(BuildPkg):
    name = "netCDF4"
    def _install(self):
        return self.sys_log("NETCDF3_DIR=%s python setup-nc3.py install --prefix=%s" % (self.blddir, self.blddir))

# FLTK (pronounced "fulltick") is a cross-platform C++ GUI toolkit.
# http://www.fltk.org/
# Version 1.1.9 required by gmsh
class Fltk(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --enable-xft --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")
    
class Python(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --enable-shared" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        self.sys_log("make install")
        BuildPkg.pypath = set_python_path(self.blddir)
    def _clean(self):
        return self.sys_log("make clean")

class Netcdf(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --with-pic -prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        # Don't use parallel make here. Breaks on some systems.
        return self.sys_log("make")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Xdmf(BuildPkg):
    name = 'Xdmf'
    def _configure(self):
        return self.sys_log("cmake %s -DCMAKE_INSTALL_PREFIX=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make")
    def _install(self):
        os.chdir(os.path.join(self.bdir, "bin"))
        return self.sys_log("install %s *.so %s" % (os.path.join(self.sdir, 'libsrc','Xdmf.py'), self.libdir))
    def _clean(self):
        return self.sys_log("make clean")

# The GNU Scientific Library (GSL) is a collection of routines for
# numerical analysis, written in C. 
# http://www.gnu.org/software/gsl/
class Gsl(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")
#    def _test(self):
#        ok = errs = 0
#        os.chdir(self.bdir)
#        logfile = self.logfile.replace('xml','txt')
#        os.system("make check > %s 2>&1" % logfile)
#        for line in open(logfile):
#            if line.find('PASS') == 0:
#                ok += 1
#            elif line.find('FAIL') == 0:
#                errs += 1
#        if errs:
#            ostr = "<Test Status=\"failed\">\n"
#        else:
#            ostr = "<Test Status=\"passed\">\n"
#        ostr += "\t<Name>gsl</Name>\n"
#        ostr += "\t<Path>%s</Path>\n" % self.sdir
#        ostr += "\t<FullName>gsl</FullName>\n"
#        ostr += "\t<FullCommandLine>make check</FullCommandLine>\n"
#        ostr += "\t<Results><Measurement><Value>"
#        ostr += cgi.escape(open(logfile).read())
#        ostr += "</Value></Measurement></Results></Test>\n"
#        f = open(self.logfile,'w')
#        f.write(ostr)
#        f.close()
#        return ok, errs

class Rlog(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --disable-docs --disable-valgrind --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Swig(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Fftw(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --enable-float --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j%s" % jobs(self.name))
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

# We don't build boost, just use the headers
class Boost(BuildPkg):
    def _install(self):
        idir = os.path.join(self.blddir, "include")
        return self.sys_log("/bin/ln -fs %s %s" % (self.bdir, idir))

class MPM(BuildPkg):
    name = "MPM"
    def _configure(self):
        e = config(self.name,'configname')
        bfile = os.path.join(self.sdir, "config", e)
        if os.path.isfile(bfile):
            os.chdir(os.path.join(self.sdir, "config"))
            self.sys_log("/bin/ln -fs %s CURRENT" % bfile)
            os.chdir(os.path.join(self.bdir, "config"))
            bfile = os.path.join(self.bdir, "config", e)            
            return self.sys_log("/bin/ln -fs %s CURRENT" % bfile)
        else:
            f = open(self.logfile, 'a')
            print >> f, "Cannot open config file %s." % bfile
            f.close()        
        return False
    def _build(self):
        # Don't use parallel make here. Breaks on some systems.
        return self.sys_log("make")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Fvm(BuildPkg):
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
        os.putenv("PYTHONPATH",os.path.join(BuildPkg.topdir, "tools","scons-local","scons-local"))
        return self.sys_log("%s/etc/buildsystem/build -j%s -C %s COMPACTOUTPUT=False PARALLEL=%s VERSION=%s COMPILER=%s" \
                                % (self.sdir, jobs(self.name), self.sdir, par, ver, comp))
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

class MEMOSA(BuildPkg):
    name = "MEMOSA"
    def _install(self):
        os.chdir(self.sdir)
        self.sys_log("install bin/* %s" % self.bindir)
        self.sys_log("install lib/* %s" % self.libdir)        
        return 0


