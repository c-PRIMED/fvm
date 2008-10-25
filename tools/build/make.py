#! /usr/bin/env python
# encoding: utf-8
#
# Memosa build script
# Martin Hunt <mmh@purdue.edu>

"""
This script configures and builds all the MEMOSA
packages and tools.

Usage: make.py [options] config[-pkgs]
Options:
  --build        Build sources. This is the default.
  --test         Run tests.
  --submit       Submit test and build results.
  --all          Build config and config-pkgs.
  -v, --verbose  Verbose output.
  --nocolor      Disable color output.

Configuration names are stored in the "config" subdirectory.
"""

import sys, cdash, time, pbs
from optparse import OptionParser
from build_packages import *
import build_utils, config, testing

def usage():
    print __doc__
    sys.exit(-1)

def main():

    parser = OptionParser()
    parser.add_option("--build (default)",
                      action="store_true", dest="build",
                      help="Build sources.")
    parser.add_option("--test",
                      action="store_true", dest="test",
                      help="Run tests.")
    parser.add_option("--submit",
                      action="store_true", dest="submit",
                      help="Submit test and build results.")
    parser.add_option("--all",
                      action="store_true", dest="all",
                      help="same as --rebuild --test --submit")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", help="verbose output")
    parser.add_option("--nocolor",
                      action="store_true", dest="nocolor", help="Disable color output.")
    parser.add_option("--nightly",
                      action="store_true", dest="nightly", help="Do nightly build and test.")

    (options, args) = parser.parse_args()
    build_utils.verbose = options.verbose

    srcpath = os.path.abspath(os.path.dirname(sys.argv[0]))

    if options.nocolor:
        build_utils.clear_colors()

    if options.nightly:
        options.build = options.test = options.submit = True

    cname = ''
    if len(args) == 1:
        cname = args[0]
    if cname == '' or not config.read(srcpath, cname, options.all):
        usage()

    if cname.endswith('-pkgs'):
        where = cname.find("-pkgs")
        cname = cname[0:where]
        
    BuildPkg.setup(cname, srcpath)
    build_utils.run_commands('before',0)
    fix_path('PATH', BuildPkg.bindir, 1, 0)
    fix_path('PYTHONPATH', BuildPkg.libdir, 1, 0)
    fix_path('LD_LIBRARY_PATH', BuildPkg.libdir, 1, 0)
    os.environ['MEMOSA_HOME'] = BuildPkg.blddir
    bs = be = ts = te = 0

    # if no options, default to build
    if not options.build and not options.test and not options.submit:
        options.build = True

    failed = 0

    if options.build:
        # first, remove all test results.  They are now invalid
        os.system("/bin/rm -f %s/*.xml" % BuildPkg.logdir)

        bs = time.time()
        open(BuildPkg.logdir+'/StartBuildTime','w').write(str(bs))
        for p in BuildPkg.packages:
            if not config.config(p.name,'Build'):
                build_utils.debug ("Skipping " + p.name)
                continue
            try:
                p.configure()
                p.build()
                p.install()
            except:
                failed = 1
                break
        be = time.time()
        open(BuildPkg.logdir+'/EndBuildTime','w').write(str(be))

    if not failed and options.test and not pbs.start(BuildPkg, cname):
        ts = time.time()
        open(BuildPkg.logdir+'/StartTestTime','w').write(str(ts))
        try:
            testing.run_all_tests(BuildPkg)
        except:
            pass
        te = time.time()
        open(BuildPkg.logdir+'/EndTestTime','w').write(str(te))

    if options.submit:
        cdash.submit(BuildPkg, cname, sys.argv, options.nightly)

    build_utils.run_commands('after',0)

    if options.build and not failed:
        f = open(os.path.join(BuildPkg.topdir, 'env.sh'), 'w')
        print >>f, "export LD_LIBRARY_PATH="+BuildPkg.libdir
        print >>f, "export PYTHONPATH="+BuildPkg.libdir
        print >>f, "export PATH=%s:$PATH" % BuildPkg.bindir
        f.close
        print "\nDONE\nYou need to do the following to use the build."
        print "export LD_LIBRARY_PATH="+BuildPkg.libdir
        print "export PYTHONPATH="+BuildPkg.libdir
        print "export PATH=%s:$PATH" % BuildPkg.bindir
        print "OR source env.sh"

    fix_path('LD_LIBRARY_PATH', BuildPkg.libdir, 1, 1)
    fix_path('PYTHONPATH', BuildPkg.libdir, 1, 1)
    fix_path('PATH', BuildPkg.bindir, 1, 1)

        
if __name__ == "__main__":
    main()
