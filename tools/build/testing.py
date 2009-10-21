#! /usr/bin/env python
# encoding: utf-8
#
# Memosa testing script
# Martin Hunt <mmh@purdue.edu>

"""
Functions for running tests
"""

import cgi
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
        if line[0] == '\n' or line[0] == '#': continue
        pdir = os.path.join(os.path.dirname(logfile), pname)
        if not os.path.isdir(pdir):
            try:
                os.makedirs(pdir)
            except:
                fatal("error creating directory " + pdir)
        tname = line.split()[0]
        if line.find('TESTDIR') >= 0:
            pdir = os.path.join(pdir, tname)
            line = line.replace('TESTDIR', pdir)
            if not os.path.isdir(pdir):
                try:
                    os.makedirs(pdir)
                except:
                    fatal("error creating directory " + pdir)

        outfile = os.path.join(pdir, tname + '.dat')
        line = line.replace('TESTOUT', outfile)
        line = line.split()
        cmd = ' '.join(line[1:])
        verbose(1, "Test %s: %s" % (tname, cmd))
        ostr = "\t<Name>%s</Name>\n" % tname
        ostr += "\t<Path>%s</Path>\n" % os.path.dirname(file)
        ostr += "\t<FullName>%s</FullName>\n" % tname
        ostr += "\t<FullCommandLine>%s</FullCommandLine>\n" % cgi.escape(cmd)
        ostr += "\t<Results>\n"
        t = time.time()
        exe = os.popen("/bin/bash -c '%s 2>&1'" % cmd)
        res = exe.read()
        ret = exe.close()
        t = time.time() - t
        ostr += "\t\t<NamedMeasurement type=\"numeric/double\" name=\"Execution Time\"><Value>%s</Value></NamedMeasurement>\n" % t
        ostr += "\t\t<Measurement><Value>"
        ostr += cgi.escape(res)
        ostr += "</Value></Measurement></Results>"
        if ret:
            ret >>= 8
            verbose_warn("%s Failed" % tname)
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
        _ok, _errs = run_tests(pname, os.path.join(root, 'TESTS'), logfile)
        ok += _ok
        errs += _errs
    return ok, errs

def run_all_tests(bld):
    errs = ok = 0

    # run any before commands
    run_commands('Testing', 'before')

    # add some dirs to our path
    tooldir = os.path.join(bld.topdir, "tools", "test")
    fix_path('PATH', tooldir, 1, 0)
    fix_path('PATH', bld.bindir, 1, 0)
    fix_path('LD_LIBRARY_PATH', bld.libdir, 1, 0)
    try:
        oldpypath = os.environ['PYTHONPATH']
    except:
        oldpypath = ''
    set_python_path(bld.blddir)

    # test each package
    for p in bld.packages:
        if not p.is_pkg:
            tok, terrs = p.test()
            ok += tok
            errs += terrs

    # restore paths
    if oldpypath:
        os.environ['PYTHONPATH'] = oldpypath
    else:
        del os.environ['PYTHONPATH']
    fix_path('LD_LIBRARY_PATH', bld.libdir, 1, 1)
    fix_path('PATH', bld.bindir, 1, 1)
    fix_path('PATH', tooldir, 1, 1)

if __name__ == "__main__":
    for testfile in find_tests(sys.argv[1]):
        print "testfile=", testfile

