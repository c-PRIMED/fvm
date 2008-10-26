#! /usr/bin/env python
# encoding: utf-8
#
# Memosa testing script
# Martin Hunt <mmh@purdue.edu>

"""
Functions for running tests
"""

import os, sys, cgi, config
from build_packages import *
from build_utils import *

# generator that finds all the test files under a directory.
def find_tests(startdir):
    for root, dirs, files in os.walk(startdir):
        for x in files:
            if x == 'TESTS':
                yield root

# run all the tests in a file. return number of errors
def run_tests(pname, file, logfile):
    errs = ok = 0
    f = open(logfile, 'a')

    for line in open(file):
        line = line.split()
        if len(line) == 0 or line[0][0] == '#': continue
        tname = line[0]
        cmd = ' '.join(line[1:])

        # special hack for mtest.  append --datafile
        if line[1].strip("\"") == "mtest.py" or line[1].strip("\"") == "lammps_cmp.py":
            pdir = os.path.join(os.path.dirname(logfile),pname)
            if not os.path.isdir(pdir):
                try:
                    os.makedirs(pdir)
                except:
                    fatal("error creating directory " + pdir)
            cmd += " --datafile %s/%s.dat" % (pdir, tname)

        debug("Test %s: %s" % (tname, cmd))
        ostr = "\t<Name>%s</Name>\n" % tname
        ostr += "\t<Path>%s</Path>\n" % os.path.dirname(file)
        ostr += "\t<FullName>%s</FullName>\n" % tname        
        ostr += "\t<FullCommandLine>%s</FullCommandLine>\n" % cgi.escape(cmd)
        ostr += "\t<Results><Measurement><Value>"
        exe = os.popen("/bin/bash -c '%s 2>&1'" % cmd)
        for line in exe:
            ostr +=  cgi.escape(line)
        ostr += "</Value></Measurement></Results>"
        ret = exe.close()
        if ret:
            ret >>= 8
            debug("%s Failed" % line[0])
            ostr = "<Test Status=\"failed\">\n" + ostr
            errs += 1
        else:
            ostr = "<Test Status=\"passed\">\n" + ostr
            ok += 1
        ostr += "</Test>\n"
        f.write(ostr)
    f.close()
    return ok, errs

# run all the tests in a directory.  return number of errors
def do_tests(pname, startdir, logfile):
    errs = ok = 0
    for root in find_tests(startdir):
        os.chdir(root)
        _ok, _errs = run_tests(pname, os.path.join(root,'TESTS'), logfile)
        ok += _ok
        errs += _errs
    return ok, errs

def run_all_tests(bp):
    errs = ok = 0

    # run any before commands
    run_commands('before', 'Testing')
    
    # add some dirs to our path
    tooldir = os.path.join(bp.topdir, "tools", "test")
    fix_path('PATH', tooldir, 1, 0)
    fix_path('PATH', bp.bindir, 1, 0)
    fix_path('LD_LIBRARY_PATH', bp.libdir, 1, 0)

    # test each package
    for p in bp.packages:
        if not config(p.name,'Build'):
            debug ("Skipping " + p.name)
            continue
        tok, terrs =  p.test()
        ok += tok
        errs += terrs

    # test integration
    tdir = os.path.join(bp.topdir,'tests')
    pmess("TEST",'MEMOSA', tdir)
    os.chdir(tdir)
    logfile = os.path.join(bp.logdir, "MEMOSA-test.xml")
    remove_file(logfile)
    tok, terrs = do_tests('MEMOSA', tdir, logfile)
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

