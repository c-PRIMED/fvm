#! /usr/bin/env python
# encoding: utf-8
#
# Memosa testing script
# Martin Hunt <mmh@purdue.edu>

"""
Functions for running tests
"""

import cgi, re
from build_packages import *
from build_utils import *
from subprocess import Popen, PIPE, STDOUT


class Test():
    def __init__(self):
        self.errs = self.ok = 0
        self.unit_errs = self.unit_ok = 0
        self.test_types = {}
        self.test_desc = {'./ptest.py': 'FVM Complex Parallel Solver Tests',
                          './testCellMark': 'TestCellMark',
                          'lammps_cmp.py': 'Full Lammps Test Runs',
                          'mpmtest.py': 'Full MPM Test Runs',
                          }

    def test_exec(self, cmd):
        #print "Executing %s" % cmd
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        res = p.stdout.read()
        rc = p.wait()

        cmd = cmd.split()[0].strip("\"")
        if cmd in self.test_types:
            self.test_types[cmd] += 1
        else:
            self.test_types[cmd] = 1

        if cmd == 'nosetests':
            ok = bad = 0
            tests = re.compile(r'Ran (.*) tests').findall(res)
            if tests:
                ok = int(tests[0])
            fails = re.compile(r'failures=(.*)\)').findall(res)
            if fails:
                bad = int(fails[0])
            self.unit_ok += ok
            self.unit_errs += bad
        else:
            if rc == 0:
                self.ok += 1
            else:
                self.errs += 1

        return rc, res

    def print_results(self):
        print '='*40
        print "PASSED  %d" % (self.ok + self.unit_ok)
        print "FAILED  %d" % (self.errs + self.unit_errs)
        print "\nTest Breakdown by type:"
        for t in self.test_types:
            if t == 'nosetests':
                print "%5d  unit tests consisting of" % self.test_types[t]
                print "       %5d PASSED" % self.unit_ok
                print "       %5d FAILED" % self.unit_errs
            else:
                if t in self.test_desc:
                    desc = self.test_desc[t]
                else:
                    desc = t
                print "%5d  %s" % (self.test_types[t], desc)

    #run all the tests in a file. return number of errors
    def run_tests(self, pname, fname, logfile):
        errs = ok = 0
        f = open(logfile, 'a')
        for line in open(fname):
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
            ostr += "\t<Path>%s</Path>\n" % os.path.dirname(fname)
            ostr += "\t<FullName>%s</FullName>\n" % tname
            ostr += "\t<FullCommandLine>%s</FullCommandLine>\n" % cgi.escape(cmd)
            ostr += "\t<Results>\n"
            t = time.time()
            err, result_text = self.test_exec(cmd)
            t = time.time() - t
            ostr += "\t\t<NamedMeasurement type=\"numeric/double\" name=\"Execution Time\"><Value>%s</Value></NamedMeasurement>\n" % t
            ostr += "\t\t<Measurement><Value>"
            ostr += cgi.escape(result_text)
            ostr += "</Value></Measurement></Results>"
            if err:
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

# generator that finds all the test files under a directory.
def find_tests(startdir):
    for root, dirs, files in os.walk(startdir):
        for x in files:
            if x == 'TESTS':
                yield root

# run all the tests in a directory.  return number of errors
def do_tests(tst, pname, startdir, logfile):
    errs = ok = 0
    for root in find_tests(startdir):
        os.chdir(root)
        _ok, _errs = tst.run_tests(pname, os.path.join(root, 'TESTS'), logfile)
        ok += _ok
        errs += _errs
    return ok, errs

def run_all_tests(bld):
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

    # Create Test Object to hold and display results
    tst = Test()

    # test each package
    for p in bld.packages:
        if not p.is_pkg:
            p.test(tst)

    tst.print_results()

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

