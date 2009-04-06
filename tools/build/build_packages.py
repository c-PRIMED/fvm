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
            Nose("pkgs/python-nose", 1),
            Mpi4py("pkgs/mpi4py", 1),
            Gsl("pkgs/gsl", 0),            
            Fltk("pkgs/fltk", 1),
            Gmsh("pkgs/gmsh", 1),
            Rlog("pkgs/rlog", 1),
            Boost("pkgs/boost", 0),
            Swig("pkgs/swig", 0),
            Fftw("pkgs/fftw", 0),
            Netcdf("pkgs/netcdf", 2),
            NetCDF4("pkgs/netCDF4", 1),
            ParMetis("pkgs/ParMetis", 1),
            Lammps("src/lammps", 1),
            MPM("src/MPM", 1),
            Fvm("src/fvm", 0),
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

class Lammps(BuildPkg):
    def _build(self):
        os.chdir('src')
        ret = self.sys_log("make -j4");
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
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

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
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")
    
class Python(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s --enable-shared" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j4")
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

# The GNU Scientific Library (GSL) is a collection of routines for
# numerical analysis, written in C. 
# http://www.gnu.org/software/gsl/
class Gsl(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j4")
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
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Swig(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Fftw(BuildPkg):
    def _configure(self):
        return self.sys_log("%s/configure --enable-float --prefix=%s" % (self.sdir, self.blddir))
    def _build(self):
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

# We don't build boost, just use the headers
class Boost(BuildPkg):
    def _install(self):
        idir = os.path.join(self.blddir, "include")
        return self.sys_log("/bin/ln -fs %s %s" % (self.sdir, idir))

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
        return self.sys_log("%s/etc/buildsystem/build -C %s COMPACTOUTPUT=False PARALLEL=%s VERSION=%s COMPILER=%s" \
                                % (self.sdir, self.sdir, par, ver, comp))
    def _install(self):
        try:
            line = os.popen("/bin/bash -c 'gcc --version 2>&1'").readline()
            vers = "gcc-" + re.compile(r'[^(]*[^)]*\) ([^\n ]*)').findall(line)[0]
        except:
            vers = '4.2.1'
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



