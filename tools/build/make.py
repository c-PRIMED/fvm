#! /usr/bin/env python
# encoding: utf-8
#
# Memosa build script
# Martin Hunt <mmh@purdue.edu>

"""
This script configures and builds all the MEMOSA
packages and tools.

Usage: make.py [options] config_name
Options:
  --clean        Clean up source directory.
  --rebuild      Remove build directory (if it exists) and build.
  --test         Run tests.
  --submit       Rubmit test and build results.
  --all          Same as --rebuild --test --submit.
  -v, --verbose  Verbose output.
  --nocolor      Disable color output.
Configuration names are stored in the "config" subdirectory.
"""

import sys
from optparse import OptionParser
from build_packages import *
import build_utils, config, testing

def usage():
    print __doc__
    sys.exit(-1)

def main(args):

    parser = OptionParser()
    parser.add_option("--build (default)",
                      action="store_true", dest="build",
                      help="Build sources.")
    parser.add_option("--rebuild",
                      action="store_true", dest="rebuild",
                      help="Remove build directory (if it exists) and build.")
    parser.add_option("--test",
                      action="store_true", dest="test",
                      help="Run tests.")
    parser.add_option("--submit",
                      action="store_true", dest="submit",
                      help="Submit test and build results.")
    parser.add_option("--all",
                      action="store_true", dest="all",
                      help="same as --rebuild --test --submit")
    parser.add_option("--clean",
                      action="store_true", dest="clean",
                      help="Clean up source directory.")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", help="verbose output")
    parser.add_option("--nocolor",
                      action="store_true", dest="nocolor", help="Disable color output.")
    (options, args) = parser.parse_args()
    build_utils.verbose = options.verbose

    if options.nocolor:
        build_utils.clear_colors()

    if options.clean:
        if len(args): usage()
        BuildPkg.setup()
        for p in BuildPkg.packages:
            p.clean()

    else:
        cname = ''
        if len(args) == 1:
            cname = args[0]
        if cname == '' or not config.read(cname):
            usage()
    
        BuildPkg.setup(cname)

        build_utils.run_commands('before',0)

        if options.all:
            options.rebuild = options.test = options.submit = True

        if options.rebuild:
            options.build = True
            os.system("/bin/rm -rf %s" % BuildPkg.blddir)
            BuildPkg.setup(cname)

        # if no options, default to build
        if not options.build and not options.test and not options.submit:
            options.build = True

        if options.build:
            for p in BuildPkg.packages:
                if config.config(p.name,'skip'):
                    build_utils.debug ("Skipping " + p.name)
                    continue
                p.configure()
                p.build()
                p.install()

        if options.test:
            testing.run_all_tests(BuildPkg)

        if options.submit:
            print "Submit not implemented yet."

        build_utils.run_commands('after',0)

        if options.build:
            print "\nDONE\nYou need to do something like the following to use the build."
            print "export LD_LIBRARY_PATH="+BuildPkg.libdir
            print "export PATH=$PATH:"+BuildPkg.bindir
        
if __name__ == "__main__":
    main(sys.argv[1:])
