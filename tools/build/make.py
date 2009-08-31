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
  -d, --debug    Debug output.
  -j num         Specify the number of jobs to run simultaneously.
  --nocolor      Disable color output.

Configuration names are stored in the "config" subdirectory.
"""

import sys, os, cdash, time, pbs, update, traceback
from optparse import OptionParser
from build_packages import BuildPkg
import build_utils, config, testing

def usage():
    print __doc__
    sys.exit(-1)

def main():

    parser = OptionParser()
    parser.add_option("--build",  action="store_true")
    parser.add_option("--test", action="store_true")
    parser.add_option("--update", action="store_true")
    parser.add_option("--submit", action="store_true")
    parser.add_option("--all", action="store_true")
    parser.add_option("-v", "--verbose", action="store_true")
    parser.add_option("-d", "--debug", action="store_true")
    parser.add_option("--nocolor", action="store_true")
    parser.add_option("--nightly", action="store_true")
    parser.add_option("--jobs", "-j")
    (options, args) = parser.parse_args()

    srcpath = os.path.abspath(os.path.dirname(sys.argv[0]))

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

    build_utils.set_options(options)

    if options.all:
        os.system("/bin/rm -rf %s" % os.path.join(os.getcwd(), "build-%s" % cname))

    BuildPkg.setup(cname, srcpath)
    build_utils.run_commands('ALL', 'before')        
    build_utils.fix_path('PATH', BuildPkg.bindir, 1, 0)
    build_utils.fix_path('LD_LIBRARY_PATH', BuildPkg.libdir, 1, 0)
    os.environ['MEMOSA_HOME'] = BuildPkg.blddir
    os.environ['MEMOSA_CONFNAME'] = cname
    build_start_time = build_end_time = test_start_time = test_end_time = 0
    try:
        oldpypath = os.environ['PYTHONPATH']
    except:
        oldpypath = ''
    BuildPkg.pypath = build_utils.set_python_path(BuildPkg.blddir)    

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
        build_start_time = time.time()
        open(BuildPkg.logdir+'/StartBuildTime','w').write(str(build_start_time))
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
                traceback.print_exc()
                build_failed = 1
                break
        build_end_time = time.time()
        open(BuildPkg.logdir+'/EndBuildTime','w').write(str(build_end_time))

        # write out env.sh
        env_name = os.path.join(os.getcwd(), 'env.sh')
        f = open(env_name, 'w')
        modules = config.config('Testing', 'modules')
        if modules:
            for m in modules.split():
                f.write('module load %s\n' % m)            
        print >> f, "export LD_LIBRARY_PATH="+BuildPkg.libdir+":$LD_LIBRARY_PATH"
        try:
            if os.environ['PYTHONPATH']:
                print >> f, "export PYTHONPATH="+os.environ['PYTHONPATH']
        except:
            pass
        print >> f, "export PATH=%s:$PATH" % BuildPkg.bindir
        print >> f, "\n# Need this to recompile MPM in its directory."
        print >> f, "export MEMOSA_CONFNAME=%s" % cname
        f.close()
        if not build_failed:
            print "\nDONE\nYou need to source %s to use this build." %  env_name)

    # TESTING
    if options.test and not pbs.start(BuildPkg, cname):
        test_start_time = time.time()
        open(BuildPkg.logdir+'/StartTestTime','w').write(str(test_start_time))
        testing.run_all_tests(BuildPkg)
        test_end_time = time.time()
        open(BuildPkg.logdir+'/EndTestTime','w').write(str(test_end_time))

    # SUBMIT
    if options.submit:
        cdash.submit(BuildPkg, cname, sys.argv, options.nightly)

    if not options.test:
        build_utils.run_commands('ALL', 'after')

    build_utils.fix_path('LD_LIBRARY_PATH', BuildPkg.libdir, 1, 1)
    if oldpypath:
        os.environ['PYTHONPATH'] = oldpypath
    else:
        del os.environ['PYTHONPATH']
    build_utils.fix_path('PATH', BuildPkg.bindir, 1, 1)

if __name__ == "__main__":
    main()
