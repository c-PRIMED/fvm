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
  --update       Update sources from the subversion repository.
  --all          Removes build directory then builds config and config-pkgs.
  -v, --verbose  Verbose output.
  --nocolor      Disable color output.

Configuration names are stored in the "config" subdirectory.
"""

import sys, cdash, time, pbs, update
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
    parser.add_option("--update",
                      action="store_true", dest="update",
                      help="Update sources from the subversion repository.")
    parser.add_option("--submit",
                      action="store_true", dest="submit",
                      help="Submit test and build results.")
    parser.add_option("--all",
                      action="store_true", dest="all",
                      help="Removes build directory then builds config and config-pkgs.")
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
        options.all = options.update = options.test = options.submit = True

    if options.all:
        options.build = True
        
    cname = ''
    if len(args) == 1:
        cname = args[0]

    packages = options.all
    if cname.endswith('-pkgs'):
        where = cname.find("-pkgs")
        cname = cname[0:where]
        packages = True

    if cname == '' or not config.read(srcpath, cname, packages):
        usage()

    if options.all:
        os.system("/bin/rm -rf %s" % os.path.join(os.getcwd(), "build-%s" % cname))

    BuildPkg.setup(cname, srcpath)
    build_utils.run_commands('ALL', 'before')        
    fix_path('PATH', BuildPkg.bindir, 1, 0)
    fix_path('LD_LIBRARY_PATH', BuildPkg.libdir, 1, 0)
    os.environ['MEMOSA_HOME'] = BuildPkg.blddir
    os.environ['MEMOSA_CONFNAME'] = cname
    bs = be = ts = te = 0
    try:
        oldpypath = os.environ['PYTHONPATH']
    except:
        oldpypath = ''
    BuildPkg.pypath = set_python_path(BuildPkg.blddir)    

    # if no options, default to build
    if not options.build and not options.test and not options.submit and not options.update:
        options.build = True

    if options.build:
        # Remove all test results.  They are now invalid
        os.system("/bin/rm -f %s/*.xml" % BuildPkg.logdir)

    # UPDATE
    if options.update:
        update.update(BuildPkg, cname, options.nightly)

    # BUILDING
    build_failed = 0
    if options.build:
        bs = time.time()
        open(BuildPkg.logdir+'/StartBuildTime','w').write(str(bs))
        for p in BuildPkg.packages:
            x = config.config(p.name,'Build')
            if x == '' or not eval(x):
                build_utils.debug ("Skipping " + p.name)
                continue
            try:
                p.configure()
                p.build()
                p.install()
            except:
                build_failed = 1
                break
        be = time.time()
        open(BuildPkg.logdir+'/EndBuildTime','w').write(str(be))

    # TESTING
    if options.test and not pbs.start(BuildPkg, cname):
        ts = time.time()
        open(BuildPkg.logdir+'/StartTestTime','w').write(str(ts))
        try:
            testing.run_all_tests(BuildPkg)
        except:
            pass
        te = time.time()
        open(BuildPkg.logdir+'/EndTestTime','w').write(str(te))

    # SUBMIT
    if options.submit:
        cdash.submit(BuildPkg, cname, sys.argv, options.nightly)

    if not options.test:
        build_utils.run_commands('ALL', 'after')

    if options.build:
        f = open(os.path.join(BuildPkg.topdir, 'env.sh'), 'w')
        print >>f, "export LD_LIBRARY_PATH="+BuildPkg.libdir+":$LD_LIBRARY_PATH"
        try:
            if os.environ['PYTHONPATH']:
                print >>f, "export PYTHONPATH="+os.environ['PYTHONPATH']
        except:
            pass
        print >>f, "export PATH=%s:$PATH" % BuildPkg.bindir
        print >>f, "\n# Need this to recompile MPM in its directory."
        print >>f, "export MEMOSA_CONFNAME=%s" % cname
        f.close
        print "\nDONE\nYou need to do the following to use the build."
        print "export LD_LIBRARY_PATH="+BuildPkg.libdir+":$LD_LIBRARY_PATH"
        try:
            if os.environ['PYTHONPATH']:
                print "export PYTHONPATH="+os.environ['PYTHONPATH']
        except:
            pass
        print "export PATH=%s:$PATH" % BuildPkg.bindir
        print "OR source env.sh"

    fix_path('LD_LIBRARY_PATH', BuildPkg.libdir, 1, 1)
    if oldpypath:
        os.environ['PYTHONPATH'] = oldpypath
    else:
        del os.environ['PYTHONPATH']
    fix_path('PATH', BuildPkg.bindir, 1, 1)

if __name__ == "__main__":
    main()
