#! /usr/bin/env python
# encoding: utf-8
#
# Memosa build script
# Martin Hunt <mmh@purdue.edu>

"""
This script configures and builds a collection of packages.

Usage: make.py [options] config
Options:
  --build        Build sources. This is the default.
  --test         Run tests. Will build first if necessary.
  --submit       Submit the latest test and build results.
  --update       Update sources from the subversion repository.
  --all          Removes build directory then builds everything.
  -v, --verbose  Verbose output.
  -d, --debug    Debug output.
  -s source_dir  Directory where the sources are located.
  -j num         Specify the number of jobs to run simultaneously.
  --nocolor      Disable color output.
  --clean        Clean up. Currently only removes old binaries from fvm directory.
  --nightly      Build, test, and submit to cdash. Flag as a "nightly" build.

Configuration names are stored in the 'config' subdirectory. That directory must be
in the source directory specified by the '-s' option (if one is given).
If there is a 'config' directory in the current directory, it will be used.
Otherwise, the build system will look for one in the same directory as 'make.py'.
"""

from build import Build
from optparse import OptionParser
from shutil import rmtree
import build_utils
import config
import testing
import sys
import os
import cdash
import time
import pbs
import moab
import update
import traceback

def usage():
    print __doc__
    sys.exit(-1)

def main():

    parser = OptionParser()
    parser.set_defaults(verbose=0)
    parser.add_option("--build", action="store_true")
    parser.add_option("--test", action="store_true")
    parser.add_option("--update", action="store_true")
    parser.add_option("--submit", action="store_true")
    parser.add_option("--all", action="store_true")
    parser.add_option("-v", "--verbose", action="count")
    parser.add_option("-d", "--debug", action="store_true")
    parser.add_option("--nocolor", action="store_true")
    parser.add_option("--nightly", action="store_true")
    parser.add_option("--clean", action="store_true")
    parser.add_option("--jobs", "-j")
    parser.add_option("-s", type="string")
    (options, args) = parser.parse_args()

    make_path = os.path.abspath(os.path.dirname(os.path.realpath(sys.argv[0])))
    cwd = os.getcwd()


    if options.nightly:
        options.update = options.test = options.submit = True

    if options.all or options.test:
        options.build = True

    cname = ''
    if len(args) == 1:
        cname = args[0]
    if cname == '':
        usage()

    sdir = ''
    if options.s:
        sdir = options.s
    else:
        if os.path.isdir(os.path.join(cwd, 'config')):
            sdir = cwd
        elif os.path.islink(sys.argv[0]):
            sdir = os.path.dirname(sys.argv[0])

    if sdir == '' or not config.read(os.path.join(sdir, 'config'), cname):
        usage()

    build_utils.set_options(options)
    if options.all:
        os.system("/bin/rm -rf %s" % os.path.join(cwd, "build-%s" % cname))

    cmd = config.config('ALL', 'before')
    if cmd and not '_MEMOSA_MAKE' in os.environ:
        cmd = ';'.join(cmd)
        os.environ['_MEMOSA_MAKE'] = '1'
        ret = os.system("/bin/bash -l -c '%s;%s'" %(cmd, ' '.join(sys.argv)))
        sys.exit(ret>>8)

    bld = Build(cname, sdir, make_path)

    # CLEAN
    if options.clean or options.all:
        for p in bld.all_packages:
            p.clean()

    build_utils.fix_path('PATH', bld.bindir, 1, 0)
    build_utils.fix_path('LD_LIBRARY_PATH', bld.libdir, 1, 0)
    os.environ['MEMOSA_HOME'] = bld.blddir
    os.environ['MEMOSA_CONFNAME'] = cname
    build_start_time = build_end_time = test_start_time = test_end_time = 0
    try:
        oldpypath = os.environ['PYTHONPATH']
    except:
        oldpypath = ''
    build_utils.set_python_path(bld.blddir)

    # if no options, default to build
    if not options.build and not options.test and not options.submit \
            and not options.update and not options.clean:
        options.build = True

    if options.build:
        # Remove all test results.  They are now invalid
        os.system("/bin/rm -f %s/*.xml" % bld.logdir)

    # UPDATE
    if options.update:
        update.update(bld, cname, options.nightly)

    # BUILDING
    build_failed = 0

    if options.build and bld.packages == []:
        print "No packages need built."

    if options.build and bld.packages != []:
        build_start_time = time.time()
        open(bld.logdir + '/StartBuildTime', 'w').write(str(build_start_time))
        for p in bld.packages:
            try:
                p.build()
            except build_utils.CompileException:
                build_failed = 1
                break
            except:
                traceback.print_exc()
                build_failed = 1
                break
        build_end_time = time.time()
        open(bld.logdir + '/EndBuildTime', 'w').write(str(build_end_time))

        # write out env.[c]sh
        env_name = build_utils.write_env(bld, cwd, cname)

        if not build_failed:
            print "\nDone with building.\nYou need to source %s to use this build.\n" % env_name

    # make sure we are back in the original directory
    os.chdir(cwd)

    # set package list for testing and submit
    bld.packages = []
    bld.build_pkg_list(False)
    bld.done()

    # TESTING
    if build_failed==0 and options.test and not pbs.start(bld, cname) and not moab.start(bld, cname):
        testdir = os.path.join(bld.logdir, 'testing')
        rmtree(testdir, True)
        os.mkdir(testdir)
        test_start_time = time.time()
        open(os.path.join(testdir, 'StartTestTime'), 'w').write(str(test_start_time))
        testing.run_all_tests(bld)
        test_end_time = time.time()
        open(os.path.join(testdir, 'EndTestTime'), 'w').write(str(test_end_time))

    # SUBMIT
    if options.submit:
        cdash.submit(bld, cname, sys.argv, options.nightly)

    if not options.test:
        build_utils.run_commands('ALL', 'after')

    build_utils.fix_path('LD_LIBRARY_PATH', bld.libdir, 1, 1)
    if oldpypath:
        os.environ['PYTHONPATH'] = oldpypath
    else:
        del os.environ['PYTHONPATH']
    build_utils.fix_path('PATH', bld.bindir, 1, 1)

    sys.exit(build_failed)

if __name__ == "__main__":
    main()
