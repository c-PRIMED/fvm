#! /usr/bin/env python
# encoding: utf-8
#
# Memosa build script
# Martin Hunt <mmh@purdue.edu>

"""
This script configures and builds all the MEMOSA
packages and tools.
"""

import sys
from optparse import OptionParser
from build_packages import *
import build_utils
import config
def main(args):


    BuildPkg.setup()

    parser = OptionParser()
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="verbose output")
    (options, args) = parser.parse_args()
    build_utils.verbose = options.verbose


    if len(args) and args[0] == "clean":
        for p in BuildPkg.packages:
            p.clean()
    elif len(args) and args[0] == "test":
        print "Tests are not implemented yet"
    else:
        if len(args) == 0 or args[0] == "" or config.read(args[0]):
            build_utils.run_commands('before',0)
            for p in BuildPkg.packages:
                if config.config(p.name,'skip'):
                    build_utils.debug ("Skipping " + p.name)
                    continue
                p.configure()
                p.build()
                p.install()
            build_utils.run_commands('after',0)
            print "\nDONE\nYou need to do something like the following to use the build."
            print "export LD_LIBRARY_PATH="+BuildPkg.libdir
            print "export PATH=$PATH:"+BuildPkg.bindir
        else:
            print "Unknown build option", args[0]
                
        
if __name__ == "__main__":
    main(sys.argv[1:])
