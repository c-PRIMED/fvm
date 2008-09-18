# Build packages
# Martin Hunt <mmh@purdue.edu>

"""
Build package definitions.
"""

import sys
from os import *
from build_utils import *
from config import config

def create_build_dir():
    for p in [BuildPkg.blddir,
              BuildPkg.logdir,
              BuildPkg.libdir,
              BuildPkg.bindir,
              ]:
        if not path.isdir(p):
            try: 
                makedirs(p)
            except:
                fatal("error creating directory " + p)

# abstract superclass
class BuildPkg:

    def setup():
        BuildPkg.topdir = getcwd()
        BuildPkg.blddir = path.join(BuildPkg.topdir, "build")
        BuildPkg.logdir = path.join(BuildPkg.blddir, "log")
        BuildPkg.bindir = path.join(BuildPkg.blddir, "bin")
        BuildPkg.libdir = path.join(BuildPkg.blddir, "lib")
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
        self.sdir = path.join(self.topdir, sdir)
        if bdir != "":
            self.bdir = path.join(self.sdir, bdir)
            if not path.isdir(self.bdir):
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
            system("tail -40 " + self.logfile)
            sys.exit(state)
        else:
            cprint('YELLOW', state)

    def configure(self):
        self.state = 'configure'
        self.logfile = path.join(self.logdir, self.name+"-conf.log")
        pmess("CONF",self.name,self.bdir)
        chdir(self.bdir)
        run_commands('before',self.name)
        self.pstatus(self._configure())

    def clean(self):
        self.state = 'clean'
        self.logfile = ''
        chdir(self.bdir)
        self._clean()
        if self.sdir != self.bdir:
            chdir(self.sdir)
            if path.isdir(self.bdir):            
                system("/bin/rm -rf %s/*" % self.bdir)

    def build(self):
        self.state = 'build'
        self.logfile = path.join(self.logdir, self.name+"-build.log")
        pmess("BUILD",self.name,self.bdir)
        chdir(self.bdir)
        self.pstatus(self._build())

    def install(self):
        self.state = 'install'
        self.logfile = path.join(self.logdir, self.name+"-install.log")
        pmess("INSTALL",self.name,self.blddir)
        chdir(self.bdir)
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
            cmd = cmd + " &>" + self.logfile
        return system("/bin/bash -c '%s'" % cmd)

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
        bindir = self.blddir+"/bin"
        return self.sys_log("install lmp_* %s" % bindir)

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
    def _build(self):
        os.putenv("PYTHONPATH",path.join(BuildPkg.topdir, "tools","scons-local","scons-local"))
        return self.sys_log("../etc/buildsystem/build -C %s" % self.sdir)

# self.__class__.__name__


