#! /usr/bin/env python
# encoding: utf-8
#
# Memosa testing script
# Martin Hunt <mmh@purdue.edu>

"""
Functions for running tests
"""

import os, sys, config
from build_packages import *
from build_utils import *

# generator that finds all the test files under a directory.
def find_tests(startdir):
    for root, dirs, files in os.walk(startdir):
        for x in files:
            if x == 'TESTS':
                yield root

# run all the tests in a file. return number of errors
def run_tests(file, logfile):
    errs = ok = 0
    for line in open(file):
        line = line.split()
        tname = line[0]
        cmd = ' '.join(line[1:])
        debug("Test %s: %s" % (tname, cmd))
        f = open(logfile, 'a')
        print >> f,"%s: EXECUTING: %s" % (tname, cmd)
        f.close()
        cmd = cmd + " >>" + logfile + " 2>&1"
        ret = os.system("/bin/bash -c '%s'" % cmd)
        if ret:
            debug("%s Failed" % line[0])
            f = open(logfile, 'a')
            print >> f, "%s: FAILED: %s" % (tname, ret)
            f.close()
            errs += 1
        else:
            ok += 1
    return ok, errs

# run all the tests in a directory.  return number of errors
def do_tests(startdir, logfile):
    errs = ok = 0
    for root in find_tests(startdir):
        os.chdir(root)
        _ok, _errs = run_tests(os.path.join(root,'TESTS'), logfile)
        ok += _ok
        errs += _errs
    return ok, errs

def run_all_tests(bp):
    errs = ok = 0

    # add some dirs to our path
    tooldir = os.path.join(bp.topdir, "tools", "test")
    fix_path('PATH', tooldir, 1, 0)
    fix_path('PATH', bp.bindir, 1, 0)
    fix_path('LD_LIBRARY_PATH', bp.libdir, 1, 0)

    # test each package
    for p in bp.packages:
        if config(p.name,'skip'):
            debug ("Skipping " + p.name)
            continue
        tok, terrs =  p.test()
        ok += tok
        errs += terrs

    # test integration
    tdir = os.path.join(bp.topdir,'tests')
    pmess("TEST",'MEMOSA', tdir)
    os.chdir(tdir)
    logfile = os.path.join(bp.logdir, "MEMOSA-test.log")
    remove_file(logfile)
    tok, terrs = do_tests(tdir, logfile)
    if terrs:
        cprint('YELLOW', "%s OK, %s FAIL" % (tok, terrs))
    else:
        cprint('GREEN', "%s OK" % tok)
    os.chdir(bp.topdir)
    errs += terrs
    ok += tok
    print '-'*40
    if errs:
        cprint('YELLOW', "%s OK, %s FAIL" % (ok, errs))
    else:
        cprint('GREEN', "%s OK" % ok)


    # restore paths
    fix_path('LD_LIBRARY_PATH', bp.libdir, 1, 1)
    fix_path('PATH', bp.bindir, 1, 1)
    fix_path('PATH', tooldir, 1, 1)

if __name__ == "__main__":
    for f in find_tests(sys.argv[1]):
        print "f=",f

