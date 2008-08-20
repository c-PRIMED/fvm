#! /usr/bin/env python
# encoding: utf-8
#
# Memosa build script
# Martin Hunt <mmh@purdue.edu>

"""
This script configures and builds all the MEMOSA
packages and tools.
"""

import sys, os
from optparse import OptionParser
from memosa_packages import *


# FIXME. Move and simplify all the directory stuff
srcdir = os.getcwd()
tmpdir = srcdir + "/build/tmp"
logdir = srcdir + "/build/log"
blddir = srcdir + "/build"
MemosaPkg.srcdir = srcdir
MemosaPkg.blddir = blddir
MemosaPkg.tmpdir = tmpdir
MemosaPkg.logdir = logdir
os.system("mkdir -p " + blddir)
os.system("mkdir -p " + logdir)
os.system("mkdir -p %s/%s" % (blddir,"bin"))
os.system("mkdir -p %s/%s" % (blddir,"lib"))
os.system("mkdir -p %s/%s" % (blddir,"lib"))
os.system("mkdir -p " + tmpdir)


parser = OptionParser()
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="verbose output")
(options, args) = parser.parse_args()


#system("tools/build/waf configure")

packages = [Fltk("pkgs/fltk"),
            Gmsh("pkgs/gmsh"),
            Lammps("src/lammps/src")]


if args == []:
    args = ["all"]


if args[0] == "clean":
    for p in packages:
        p.clean()
elif args[0] == "test":
    print "Tests are not implemented yet"
elif args[0] == "all":
    for p in packages:
        p.configure()
        p.build()
        p.install()
else:
    print "Unknown build option", args[0]

# FIXME
if tmpdir != '':
    os.system("rm -rf " + tmpdir)

