# Memosa packages
# Martin Hunt <mmh@purdue.edu>

"""
Memosa package definitions.
"""

import sys
from os import *
from build_utils import *

# abstract superclass
class MemosaPkg:
    def __init__(self, dir):
        self.dir = self.srcdir+'/'+dir

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
        self.logfile = self.logdir+"/"+self.name+"-conf.log"
        pmess("CONF",self.name,self.dir)
        chdir(self.dir)
        self.pstatus(self._configure())
        chdir(self.blddir);

    def clean(self):
        chdir(self.dir)
        self._clean()
        chdir(self.blddir);

    def build(self):
        self.logfile = self.logdir+"/"+self.name+"-build.log"
        pmess("BUILD",self.name,self.dir)
        chdir(self.dir)
        self.pstatus(self._build())
        chdir(self.blddir);

    def install(self):
        self.logfile = self.logdir+"/"+self.name+"-install.log"
        pmess("INSTALL",self.name,self.blddir)
        chdir(self.dir)
        self.pstatus(self._install())
        chdir(self.blddir);

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
# Packages
#########################################################

class Lammps(MemosaPkg):
    name = "lammps"
    def _build(self):
        return system("make openmpi &>" + self.logfile)
    def _clean(self):
        system("make clean-all")
        system("/bin/rm -f lmp_*")
        return 0
    def _install(self):
        bindir = self.blddir+"/bin"
        return system("install lmp_* %s &>%s" % (bindir,self.logfile))

class Gmsh(MemosaPkg):
    name = "gmsh"
    def _configure(self):
        return system("./configure --prefix=%s &>%s" % (self.blddir,self.logfile))
    def _build(self):
        return system("make &>" + self.logfile)
    def _install(self):
        return system("make install &>" + self.logfile)
    def _clean(self):
        return system("make clean")
        
    



