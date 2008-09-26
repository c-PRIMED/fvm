# Build packages
# Martin Hunt <mmh@purdue.edu>

"""
Build package definitions.
"""

import sys, os
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

    def setup(cname=''):
        BuildPkg.topdir = os.getcwd()
        if cname:
            cname = '-' + cname
        BuildPkg.blddir = os.path.join(BuildPkg.topdir, "build%s" % cname)
        BuildPkg.logdir = os.path.join(BuildPkg.blddir, "log")
        BuildPkg.bindir = os.path.join(BuildPkg.blddir, "bin")
        BuildPkg.libdir = os.path.join(BuildPkg.blddir, "lib")
        create_build_dir()

        BuildPkg.packages = [Gsl("pkgs/gsl", "build"),
                Fltk("pkgs/fltk"),
                Gmsh("pkgs/gmsh"),
                Rlog("pkgs/rlog"),
                Lammps("src/lammps/src"),
                Fvm("src/fvm","build"),
                ]

    setup = staticmethod(setup)

    def __init__(self, sdir, bdir=""):
        self.sdir = os.path.join(self.topdir, sdir)
        if bdir != "":
            self.bdir = os.path.join(self.sdir, bdir)
            if not os.path.isdir(self.bdir):
                try: 
                    makedirs(self.bdir)
                except:
                    fatal("error creating directory " + self.bdir)
        else:
            self.bdir = self.sdir

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
        os.chdir(self.bdir)
        run_commands('before',self.name)
        self.pstatus(self._configure())

    def clean(self):
        self.state = 'clean'
        self.logfile = ''
        os.chdir(self.bdir)
        self._clean()
        if self.sdir != self.bdir:
            os.chdir(self.sdir)
            if path.isdir(self.bdir):            
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
        run_commands('after',self.name)

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

#########################################################
# PACKAGES 
#
# These contain the build methods for each package.
# These will mostly be the same, but we implement
# them this way to allow maximum flexibility.
#########################################################

class Lammps(BuildPkg):
    name = "lammps"
    def _build(self):
        return self.sys_log("make -j4");
    def _clean(self):
        self.sys_log("make clean-all")
        self.sys_log("/bin/rm -f lmp_*")
        return 0
    def _install(self):
        self.sys_log("install lmp_* %s" % self.bindir)
        # for testing, we want a consistent name for the executable
        return self.sys_log("/bin/ln -fs %s/lmp_* %s/lammps" % (self.bindir, self.bindir))

class Gmsh(BuildPkg):
    name = "gmsh"
    def _configure(self):
        return self.sys_log("./configure --prefix=%s --with-gsl-prefix=%s" % (self.blddir, self.blddir))
    def _build(self):
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

# FLTK (pronounced "fulltick") is a cross-platform C++ GUI toolkit.
# http://www.fltk.org/
# Version 1.1.9 required by gmsh
class Fltk(BuildPkg):
    name = "fltk"
    def _configure(self):
        return self.sys_log("./configure --prefix=%s" % self.blddir)
    def _build(self):
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")
        
# The GNU Scientific Library (GSL) is a collection of routines for
# numerical analysis, written in C. 
# http://www.gnu.org/software/gsl/
class Gsl(BuildPkg):
    name = "gsl"
    def _configure(self):
        return self.sys_log("../configure --prefix=%s" % self.blddir)
    def _build(self):
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")
    

class Rlog(BuildPkg):
    name = "rlog"
    def _configure(self):
        return self.sys_log("./configure --disable-docs --disable-valgrind  --prefix=%s" % self.blddir)
    def _build(self):
        return self.sys_log("make -j4")
    def _install(self):
        return self.sys_log("make install")
    def _clean(self):
        return self.sys_log("make clean")

class Fvm(BuildPkg):
    name = "fvm"
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
        return self.sys_log("/bin/ln -fs %s %s" % (self.blddir, pdir))
    def _build(self):
        os.putenv("PYTHONPATH",os.path.join(BuildPkg.topdir, "tools","scons-local","scons-local"))
        return self.sys_log("../etc/buildsystem/build -C %s" % self.sdir)
    def _install(self):
        try:
            line = os.popen("/bin/bash -c 'gcc --version 2>&1'").readline()
            vers = "gcc-" + line.split()[2]
        except:
            vers = '4.2.1'
        pdir = os.path.join(self.sdir, "build", self.getArch(), vers, "debug", "bin")
        os.chdir(pdir)
        self.sys_log("install testLinearSolver %s" % self.bindir)
        return self.sys_log("install -t %s *.so" % self.libdir)


# self.__class__.__name__


