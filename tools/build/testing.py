#! /usr/bin/env python
# encoding: utf-8
#
# Memosa testing script
# Martin Hunt <mmh@purdue.edu>

"""
Functions for running tests
"""

import re
from build_packages import *
from build_utils import *
from subprocess import Popen, PIPE, STDOUT
from xml.sax.saxutils import escape
import xml.dom.minidom

class JTest:
    def __init__(self):
        self.ok = 0
        self.errs = 0

    def print_results(self):
        print '='*40
        print "PASSED  %d" % (self.ok)
        print "FAILED  %d" % (self.errs)

    def test_exec(self, cmd):
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT)
        res = p.stdout.read()
        rc = p.wait()
        return rc, res

    def fix_nose_xml(self, fname, pname):
        dom = xml.dom.minidom.parse(fname)
        ts = dom.getElementsByTagName("testsuite")[0]
        cases = dom.getElementsByTagName("testcase")
        for case in cases:
            cls = case.getAttribute("classname")
            case.setAttribute("classname", "%s.%s" % (pname,cls))
        fp = open(fname, 'w')
        fp.write(dom.toxml())
        fp.close

    def get_dom(self, pname, logdir):
        xmlname = os.path.join(logdir, "%s-test.xml" % pname)
        try:
            dom = xml.dom.minidom.parse(xmlname)
        except IOError:
            dom = xml.dom.minidom.parseString('<testsuite />')
        return dom

    #run all the tests in a file. return number of errors
    def run_tests(self, pname, fname, logdir):
        if not os.path.isdir(logdir):
            try:
                os.makedirs(logdir)
            except:
                fatal("error creating directory " + logdir)
        errs = ok = 0
        for line in open(fname):
            if line[0] == '\n' or line[0] == '#': continue
            tname = line.split()[0]
            pdir=''
            if line.find('TESTDIR') >= 0 or line.find('TESTOUT') >= 0:
                pdir = os.path.join(logdir, tname)
                if not os.path.isdir(pdir):
                    try:
                        os.makedirs(pdir)
                    except:
                        fatal("error creating directory " + pdir)
            line = line.replace('TESTDIR', pdir)
            line = line.replace('TESTOUT', os.path.join(pdir, tname + '.dat'))
            line = line.split()

            # force nosetests --with-xunit
            if line[1] == 'nosetests':
                nose = True
                xfile = os.path.join(logdir, '%s-test.xml' % tname)
                line[1] = 'nosetests --nologcapture --with-xunit --xunit-file=%s' % xfile
            else:
                nose = False
            cmd = ' '.join(line[1:])
            verbose(1, "Test %s: %s" % (tname, cmd))
            t = time.time()
            err, result_text = self.test_exec(cmd)
            t = time.time() - t
            if err:
                verbose_warn("%s Failed" % tname)
                errs += 1
                self.errs += 1
            else:
                ok += 1
                self.ok += 1
            if nose:
                self.fix_nose_xml(xfile, pname)
            else:
                # add testcase to the dom
                ts = self.dom.getElementsByTagName("testsuite")[0]
                tc = self.dom.createElement('testcase')
                tc.setAttribute('classname', pname)
                #tc.setAttribute('name', '%s.%s' % (pname,tname))
                tc.setAttribute('name', '%s.%s' % (pname,tname))
                tc.setAttribute('time', str(t))
                if err:
                    f = self.dom.createElement('failure')
                    f.setAttribute('type', 'unknown')
                    f.setAttribute('message', result_text)
                    tc.appendChild(f)
                ts.appendChild(tc)
        return ok, errs

class Test:
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
            else:
                fails = re.compile(r'errors=(.*)\)').findall(res)
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
        logfile = logfile + '-test.xml'
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
            ostr += "\t<FullCommandLine>%s</FullCommandLine>\n" % escape(cmd)
            ostr += "\t<Results>\n"
            t = time.time()
            err, result_text = self.test_exec(cmd)
            t = time.time() - t
            ostr += "\t\t<NamedMeasurement type=\"numeric/double\" name=\"Execution Time\"><Value>%s</Value></NamedMeasurement>\n" % t
            ostr += "\t\t<Measurement><Value>"
            # nosetests puts in an escape sequence we need to strip
            # before xml quoting.
            result_text = escape(result_text.rstrip("\x1b[?1034h"))
            ostr += result_text
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
def do_tests(tst, pname, startdir, logdir):
    if isinstance(tst, JTest):
        tst.dom = tst.get_dom(pname, logdir)
    errs = ok = 0

    for root in find_tests(startdir):
        os.chdir(root)
        _ok, _errs = tst.run_tests(pname, os.path.join(root, 'TESTS'), logdir)
        ok += _ok
        errs += _errs
    if isinstance(tst, JTest) and (ok or errs):
        ts = tst.dom.getElementsByTagName("testsuite")[0]
        ts.setAttribute('tests', str(ok+errs))
        ts.setAttribute('errors', '0')
        ts.setAttribute('failures', str(errs))
        ts.setAttribute('skip', '0')
        xmlname = os.path.join(logdir, "%s-test.xml" % pname)
        fp = open(xmlname, 'w')
        fp.write(tst.dom.toxml())
        fp.close
    return ok, errs

def run_all_tests(bld, ttype):
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
    if ttype == 'dash':
        tst = Test()
    else:
        tst = JTest()
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

